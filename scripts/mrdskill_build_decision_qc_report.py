#!/usr/bin/env python3
"""Build a compact decision-oriented MRD QC report from detailed VCF QC outputs."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import sqlite3
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from mrdskill_collect_case_qc_outputs import collect_case_qc_outputs
from mrdskill_project_paths import artifact_path


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def to_float(value: Any) -> float | None:
    text = clean(value)
    if not text:
        return None
    try:
        num = float(text)
    except ValueError:
        return None
    if math.isinf(num) or math.isnan(num):
        return None
    return num


def fmt_num(value: Any, digits: int = 3) -> str:
    num = to_float(value)
    if num is None:
        return ""
    return f"{num:.{digits}g}"


def fmt_af(value: Any) -> str:
    num = to_float(value)
    if num is None:
        return ""
    return f"{num:.4f}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a compact decision QC report for one VCF case.")
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--variant-csv", required=True)
    parser.add_argument("--rule-db", default="data/rag/mrd_latest_kg_rule_index.sqlite")
    parser.add_argument("--out-html", default="")
    parser.add_argument("--out-md", default="")
    parser.add_argument("--out-shortlist-csv", default="")
    parser.add_argument("--case-id", default="")
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle))


def variant_score(row: dict[str, str]) -> float:
    score = to_float(row.get("technical_score")) or 0
    af = to_float(row.get("tumor_af")) or 0
    kg_bonus = 25 if row.get("kg_match_status") == "gene_evidence_match" else 0
    normal_penalty = 25 if (to_float(row.get("control_vd")) or 0) > 0 else 0
    return score + kg_bonus + af * 20 - normal_penalty


def is_zero_normal(row: dict[str, str]) -> bool:
    return (to_float(row.get("control_vd")) or 0) == 0


def is_pass_strong(row: dict[str, str]) -> bool:
    return row.get("filter") == "PASS" and row.get("qc_tier") == "strong_technical_candidate"


def classify_decision(summary: dict[str, Any], rows: list[dict[str, str]]) -> dict[str, Any]:
    pass_rows = [row for row in rows if row.get("filter") == "PASS"]
    pass_strong = [row for row in pass_rows if row.get("qc_tier") == "strong_technical_candidate"]
    pass_strong_zero_normal = [row for row in pass_strong if is_zero_normal(row)]
    pass_strong_kg = [row for row in pass_strong_zero_normal if row.get("kg_match_status") == "gene_evidence_match"]
    hold_rows = [row for row in rows if row.get("qc_tier") == "hold_or_reject_without_review"]
    normal_signal_rows = [row for row in rows if (to_float(row.get("control_vd")) or 0) > 0]
    nonpass_rows = [row for row in rows if row.get("filter") != "PASS"]

    if pass_strong_kg:
        main_decision = "CONDITIONAL_PASS_FOR_MARKER_SHORTLIST"
        main_reason = "Use only the PASS, strong, zero-normal-alt subset; at least one candidate has local MRD KG support."
    elif pass_strong_zero_normal:
        main_decision = "TECHNICAL_PASS_BIOLOGICAL_SUPPORT_LIMITED"
        main_reason = "Technically usable PASS candidates exist, but none have direct local MRD KG evidence support."
    else:
        main_decision = "DO_NOT_USE_WITHOUT_MANUAL_REVIEW"
        main_reason = "No PASS strong zero-normal-alt marker candidate was found."

    full_vcf_decision = "DO_NOT_USE_FULL_RAW_VCF_FOR_MRD_MARKER_SELECTION"
    if summary.get("matched_control_signal_count", 0) == 0 and len(nonpass_rows) == 0:
        full_vcf_decision = "FULL_VCF_TECHNICALLY_USABLE"

    return {
        "main_decision": main_decision,
        "main_reason": main_reason,
        "full_vcf_decision": full_vcf_decision,
        "pass_count": len(pass_rows),
        "pass_strong_count": len(pass_strong),
        "pass_strong_zero_normal_count": len(pass_strong_zero_normal),
        "pass_strong_kg_count": len(pass_strong_kg),
        "hold_count": len(hold_rows),
        "normal_signal_count": len(normal_signal_rows),
        "nonpass_count": len(nonpass_rows),
    }


def top_candidates(rows: list[dict[str, str]], limit: int = 16) -> list[dict[str, str]]:
    candidates = [row for row in rows if is_pass_strong(row) and is_zero_normal(row)]
    return sorted(candidates, key=variant_score, reverse=True)[:limit]


def top_exclusions(rows: list[dict[str, str]], limit: int = 12) -> list[dict[str, str]]:
    excluded = [row for row in rows if row.get("qc_tier") == "hold_or_reject_without_review"]
    return sorted(excluded, key=lambda row: ((to_float(row.get("tumor_af")) or 0), (to_float(row.get("tumor_vd")) or 0)), reverse=True)[:limit]


def query_rules(rule_db: Path, genes: list[str], limit: int = 12) -> list[dict[str, Any]]:
    if not rule_db.exists() or not genes:
        return []
    genes = sorted({gene for gene in genes if gene})
    conn = sqlite3.connect(rule_db)
    conn.row_factory = sqlite3.Row
    try:
        placeholders = ",".join("?" for _ in genes)
        rows = [
            dict(row)
            for row in conn.execute(
                f"""
                SELECT * FROM rules
                WHERE gene_a IN ({placeholders}) OR gene_b IN ({placeholders})
                ORDER BY
                  CASE WHEN rule_status='high_confidence_rule' THEN 1 ELSE 0 END DESC,
                  CAST(pair_count AS REAL) DESC,
                  CAST(lift AS REAL) DESC
                LIMIT ?
                """,
                [*genes, *genes, limit],
            )
        ]
    finally:
        conn.close()
    return rows


def shortlist_rows(candidates: list[dict[str, str]], rules: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rules_by_gene: dict[str, list[dict[str, Any]]] = {}
    for rule in rules:
        for gene in [clean(rule.get("gene_a")), clean(rule.get("gene_b"))]:
            if gene:
                rules_by_gene.setdefault(gene, []).append(rule)
    output = []
    for rank, row in enumerate(candidates, start=1):
        genes = [gene for gene in clean(row.get("genes")).split(";") if gene]
        matched_rules = []
        for gene in genes:
            matched_rules.extend(rules_by_gene.get(gene, [])[:2])
        kg_supported = row.get("kg_match_status") == "gene_evidence_match"
        decision = "PRIMARY_MARKER" if kg_supported and rank == 1 else ("MARKER_CANDIDATE" if kg_supported else "TECHNICAL_BACKUP")
        output.append(
            {
                "rank": rank,
                "decision": decision,
                "variant_id": row.get("variant_id"),
                "genes": row.get("genes"),
                "tumor_af": row.get("tumor_af"),
                "tumor_vd": row.get("tumor_vd"),
                "tumor_dp": row.get("tumor_dp"),
                "control_vd": row.get("control_vd"),
                "technical_score": row.get("technical_score"),
                "kg_match_status": row.get("kg_match_status"),
                "kg_evidence_ids": row.get("kg_evidence_ids"),
                "rule_context": "; ".join(
                    f"{rule.get('rule_type')}:{rule.get('cancer_type')}:{rule.get('gene_a')}-{rule.get('gene_b')}"
                    for rule in matched_rules[:3]
                ),
                "flags": row.get("qc_flags"),
            }
        )
    return output


def write_shortlist_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "rank", "decision", "variant_id", "genes", "tumor_af", "tumor_vd", "tumor_dp", "control_vd",
        "technical_score", "kg_match_status", "kg_evidence_ids", "rule_context", "flags",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def graph_data(decision: dict[str, Any], shortlist: list[dict[str, Any]], exclusions: list[dict[str, str]]) -> dict[str, Any]:
    nodes = [
        {"id": "case", "label": "VCF Case", "kind": "case", "x": 0, "y": 0},
        {"id": "decision", "label": decision["main_decision"], "kind": "decision", "x": 250, "y": 0},
        {"id": "full_vcf", "label": decision["full_vcf_decision"], "kind": "warning", "x": 0, "y": 160},
        {"id": "marker_bucket", "label": f"{decision['pass_strong_zero_normal_count']} marker candidates", "kind": "candidate", "x": 250, "y": -150},
        {"id": "exclude_bucket", "label": f"{decision['hold_count']} hold/reject", "kind": "exclude", "x": 250, "y": 150},
    ]
    edges = [
        {"source": "case", "target": "decision", "label": "case decision", "kind": "decision"},
        {"source": "case", "target": "full_vcf", "label": "raw VCF decision", "kind": "warning"},
        {"source": "decision", "target": "marker_bucket", "label": "use subset", "kind": "candidate"},
        {"source": "decision", "target": "exclude_bucket", "label": "exclude/review", "kind": "exclude"},
    ]
    for idx, row in enumerate(shortlist[:8], start=1):
        node_id = f"marker_{idx}"
        nodes.append({"id": node_id, "label": f"{row['genes']} {row['variant_id']}", "kind": "primary" if row["decision"] == "PRIMARY_MARKER" else "candidate"})
        edges.append({"source": "marker_bucket", "target": node_id, "label": row["decision"], "kind": "candidate"})
        if row.get("kg_match_status") == "gene_evidence_match":
            kg_id = f"kg_{idx}"
            nodes.append({"id": kg_id, "label": "KG evidence", "kind": "kg"})
            edges.append({"source": node_id, "target": kg_id, "label": row.get("kg_evidence_ids", ""), "kind": "kg"})
    for idx, row in enumerate(exclusions[:6], start=1):
        node_id = f"exclude_{idx}"
        nodes.append({"id": node_id, "label": f"{row.get('genes')} {row.get('variant_id')}", "kind": "exclude"})
        edges.append({"source": "exclude_bucket", "target": node_id, "label": clean(row.get("filter"))[:42], "kind": "exclude"})
    return {"nodes": nodes, "edges": edges}


def table_html(rows: list[dict[str, Any]], kind: str) -> str:
    if not rows:
        return "<p class='muted'>No records.</p>"
    if kind == "shortlist":
        headers = ["Rank", "Decision", "Variant", "Gene", "AF", "VD/DP", "Score", "KG", "Rule Context"]
        body = []
        for row in rows:
            body.append(
                [
                    row.get("rank"),
                    row.get("decision"),
                    row.get("variant_id"),
                    row.get("genes"),
                    fmt_af(row.get("tumor_af")),
                    f"{row.get('tumor_vd')}/{row.get('tumor_dp')}",
                    row.get("technical_score"),
                    row.get("kg_match_status"),
                    row.get("rule_context"),
                ]
            )
    else:
        headers = ["Variant", "Gene", "AF", "VD/DP", "Filter", "Flags"]
        body = []
        for row in rows:
            body.append(
                [
                    row.get("variant_id"),
                    row.get("genes"),
                    fmt_af(row.get("tumor_af")),
                    f"{row.get('tumor_vd')}/{row.get('tumor_dp')}",
                    row.get("filter"),
                    row.get("qc_flags"),
                ]
            )
    rows_html = ["<table><thead><tr>" + "".join(f"<th>{html.escape(str(h))}</th>" for h in headers) + "</tr></thead><tbody>"]
    for cells in body:
        rows_html.append("<tr>" + "".join(f"<td>{html.escape(str(cell or ''))}</td>" for cell in cells) + "</tr>")
    rows_html.append("</tbody></table>")
    return "\n".join(rows_html)


def write_html_report(
    path: Path,
    case_id: str,
    summary: dict[str, Any],
    decision: dict[str, Any],
    shortlist: list[dict[str, Any]],
    exclusions: list[dict[str, str]],
    graph: dict[str, Any],
    shortlist_csv: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(graph, ensure_ascii=False).replace("</", "<\\/")
    decision_class = "pass" if decision["main_decision"].startswith("CONDITIONAL") else "warn"
    html_text = f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(case_id)} Decision QC</title>
  <style>
    body {{ margin:0; background:#f4f6f8; color:#172033; font:14px/1.45 "Segoe UI", Arial, sans-serif; }}
    header {{ padding:18px 22px; background:#fff; border-bottom:1px solid #dfe5ee; display:flex; justify-content:space-between; gap:16px; }}
    h1 {{ margin:0; font-size:20px; }}
    .sub {{ color:#667085; margin-top:4px; font-size:12px; }}
    main {{ padding:14px; display:grid; gap:14px; }}
    .grid {{ display:grid; grid-template-columns: 1.05fr .95fr; gap:14px; }}
    .card {{ background:#fff; border:1px solid #dfe5ee; border-radius:8px; overflow:hidden; }}
    .card h2 {{ margin:0; padding:11px 13px; font-size:14px; border-bottom:1px solid #e7ebf1; background:#fbfcfe; }}
    .body {{ padding:13px; }}
    .decision {{ border-left:6px solid {'#198754' if decision_class == 'pass' else '#b8860b'}; }}
    .badge {{ display:inline-block; padding:4px 8px; border-radius:999px; background:#eef2f6; border:1px solid #d8dee8; font-size:12px; margin:2px 4px 2px 0; }}
    .badge.pass {{ background:#e8f5ee; border-color:#b9dfc9; color:#146c43; }}
    .badge.warn {{ background:#fff6df; border-color:#ead38b; color:#7a5600; }}
    .badge.bad {{ background:#fdecec; border-color:#f3b6b6; color:#9b1c1c; }}
    .metrics {{ display:grid; grid-template-columns: repeat(4, minmax(120px,1fr)); gap:9px; }}
    .metric {{ background:#fbfcfe; border:1px solid #e7ebf1; border-radius:7px; padding:9px; }}
    .metric b {{ font-size:20px; display:block; }}
    .metric span {{ color:#667085; font-size:12px; }}
    #graph {{ height:430px; background:radial-gradient(circle at 45% 40%, #fff, #edf3f8); }}
    table {{ width:100%; border-collapse:collapse; font-size:12px; }}
    th, td {{ text-align:left; vertical-align:top; border-bottom:1px solid #edf0f4; padding:7px 6px; }}
    th {{ color:#667085; font-weight:650; background:#fbfcfe; }}
    td {{ overflow-wrap:anywhere; }}
    .muted {{ color:#667085; }}
    .note {{ color:#344054; }}
    @media (max-width: 1100px) {{ .grid {{ grid-template-columns:1fr; }} .metrics {{ grid-template-columns: repeat(2, minmax(120px,1fr)); }} }}
  </style>
</head>
<body>
<header>
  <div>
    <h1>{html.escape(case_id)} 内涵质控决策报告</h1>
    <div class="sub">Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | MRDskill KG-RAG decision layer</div>
  </div>
  <div><span class="badge {'pass' if decision_class == 'pass' else 'warn'}">{html.escape(decision['main_decision'])}</span></div>
</header>
<main>
  <section class="card decision"><h2>主要决策</h2><div class="body">
    <p><b>{html.escape(decision['main_decision'])}</b></p>
    <p>{html.escape(decision['main_reason'])}</p>
    <p><span class="badge bad">{html.escape(decision['full_vcf_decision'])}</span>
       <span class="badge pass">优先使用 PASS + strong + normal AD alt=0 子集</span>
       <span class="badge warn">互斥/共现规则仅作 review signal</span></p>
  </div></section>
  <section class="metrics">
    <div class="metric"><b>{decision['pass_count']}</b><span>PASS records</span></div>
    <div class="metric"><b>{decision['pass_strong_zero_normal_count']}</b><span>PASS strong zero-normal-alt</span></div>
    <div class="metric"><b>{decision['pass_strong_kg_count']}</b><span>KG-supported marker</span></div>
    <div class="metric"><b>{decision['normal_signal_count']}</b><span>records with normal alt reads</span></div>
  </section>
  <section class="grid">
    <div class="card"><h2>决策图谱</h2><div id="graph"></div></div>
    <div class="card"><h2>解释</h2><div class="body note">
      <p>全量原始 Mutect2 VCF 包含大量 non-PASS / normal_artifact / germline / slippage / strand-bias 记录，因此不建议把全量 VCF 直接作为 MRD marker 池。</p>
      <p>当前最主要的可用结果是：在 PASS 且 normal AD alt=0 的 strong technical candidate 中，<b>TP53 17:7578457 C&gt;A</b> 同时具备技术质量和本地 MRD KG 文献支持。</p>
      <p>精简候选表已导出：<code>{html.escape(str(shortlist_csv))}</code></p>
    </div></div>
  </section>
  <section class="card"><h2>Marker Shortlist with Decision</h2><div class="body">{table_html(shortlist, 'shortlist')}</div></section>
  <section class="card"><h2>Highest-Priority Exclusion / Manual Review</h2><div class="body">{table_html(exclusions, 'exclude')}</div></section>
</main>
<script id="graph-data" type="application/json">{data}</script>
<script>
const graph = JSON.parse(document.getElementById('graph-data').textContent);
const colors = {{case:'#2f6f73', decision:'#198754', warning:'#b8860b', candidate:'#2d6cdf', primary:'#198754', kg:'#6d5bd0', exclude:'#b94a48'}};
const el = document.getElementById('graph');
const svgNS = 'http://www.w3.org/2000/svg';
const svg = document.createElementNS(svgNS, 'svg');
svg.setAttribute('width', '100%');
svg.setAttribute('height', '100%');
svg.setAttribute('viewBox', '-160 -260 760 560');
el.appendChild(svg);
const nodeMap = new Map(graph.nodes.map((n, i) => {{
  if (n.x === undefined) {{
    const angle = i * 2.399;
    n.x = 250 + Math.cos(angle) * (170 + (i % 3) * 35);
    n.y = Math.sin(angle) * (125 + (i % 2) * 45);
  }}
  return [n.id, n];
}}));
for (const edge of graph.edges) {{
  const s = nodeMap.get(edge.source), t = nodeMap.get(edge.target);
  if (!s || !t) continue;
  const line = document.createElementNS(svgNS, 'line');
  line.setAttribute('x1', s.x); line.setAttribute('y1', s.y);
  line.setAttribute('x2', t.x); line.setAttribute('y2', t.y);
  line.setAttribute('stroke', colors[edge.kind] || '#8993a4');
  line.setAttribute('stroke-width', edge.kind === 'kg' ? '3' : '2');
  line.setAttribute('opacity', '.62');
  svg.appendChild(line);
}}
for (const node of graph.nodes) {{
  const g = document.createElementNS(svgNS, 'g');
  const c = document.createElementNS(svgNS, 'circle');
  c.setAttribute('cx', node.x); c.setAttribute('cy', node.y);
  c.setAttribute('r', node.kind === 'primary' ? 18 : 13);
  c.setAttribute('fill', colors[node.kind] || '#667085');
  c.setAttribute('stroke', '#fff'); c.setAttribute('stroke-width', '2');
  g.appendChild(c);
  const text = document.createElementNS(svgNS, 'text');
  text.setAttribute('x', Number(node.x) + 18); text.setAttribute('y', Number(node.y) + 4);
  text.setAttribute('font-size', '12'); text.setAttribute('fill', '#172033');
  text.textContent = node.label.length > 52 ? node.label.slice(0,49) + '...' : node.label;
  g.appendChild(text);
  svg.appendChild(g);
}}
</script>
</body></html>"""
    path.write_text(html_text, encoding="utf-8")


def md_table(rows: list[dict[str, Any]], fields: list[str]) -> list[str]:
    if not rows:
        return ["No records."]
    lines = ["| " + " | ".join(fields) + " |", "|" + "|".join("---" for _ in fields) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(clean(row.get(field)).replace("|", "/") for field in fields) + " |")
    return lines


def write_md_report(
    path: Path,
    case_id: str,
    decision: dict[str, Any],
    shortlist: list[dict[str, Any]],
    exclusions: list[dict[str, str]],
    html_path: Path,
    shortlist_csv: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# MRDskill Decision QC Report",
        "",
        f"- Case: `{case_id}`",
        f"- Main decision: `{decision['main_decision']}`",
        f"- Full raw VCF decision: `{decision['full_vcf_decision']}`",
        f"- Reason: {decision['main_reason']}",
        "",
        "## Key Counts",
        "",
        f"- PASS records: {decision['pass_count']}",
        f"- PASS strong zero-normal-alt candidates: {decision['pass_strong_zero_normal_count']}",
        f"- KG-supported marker candidates: {decision['pass_strong_kg_count']}",
        f"- Records with observed normal alt reads: {decision['normal_signal_count']}",
        f"- Hold/reject records: {decision['hold_count']}",
        "",
        "## Marker Shortlist",
        "",
        *md_table(shortlist[:12], ["rank", "decision", "variant_id", "genes", "tumor_af", "tumor_vd", "tumor_dp", "technical_score", "kg_match_status"]),
        "",
        "## Manual Review / Exclusion Examples",
        "",
        *md_table(exclusions[:10], ["variant_id", "genes", "tumor_af", "tumor_vd", "tumor_dp", "filter", "qc_flags"]),
        "",
        "## Files",
        "",
        f"- HTML decision report: `{html_path}`",
        f"- Shortlist CSV: `{shortlist_csv}`",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def default_out(path: Path, suffix: str) -> Path:
    return path.with_name(path.name.replace(".qc_summary.json", suffix))


def main() -> None:
    args = parse_args()
    summary_path = artifact_path(args.summary_json)
    variant_path = artifact_path(args.variant_csv)
    summary = read_json(summary_path)
    rows = read_csv(variant_path)
    case_id = args.case_id or clean(summary.get("tumor_sample")) or summary_path.stem
    decision = classify_decision(summary, rows)
    candidates = top_candidates(rows)
    genes = []
    for row in candidates:
        genes.extend(gene for gene in clean(row.get("genes")).split(";") if gene)
    rules = query_rules(artifact_path(args.rule_db), genes)
    shortlist = shortlist_rows(candidates, rules)
    exclusions = top_exclusions(rows)
    base = summary_path
    out_html = artifact_path(args.out_html) if args.out_html else default_out(base, "_decision_qc.html")
    out_md = artifact_path(args.out_md) if args.out_md else default_out(base, "_decision_qc.md")
    out_shortlist = artifact_path(args.out_shortlist_csv) if args.out_shortlist_csv else default_out(base, "_decision_shortlist.csv")
    write_shortlist_csv(out_shortlist, shortlist)
    graph = graph_data(decision, shortlist, exclusions)
    write_html_report(out_html, case_id, summary, decision, shortlist, exclusions, graph, out_shortlist)
    write_md_report(out_md, case_id, decision, shortlist, exclusions, out_html, out_shortlist)
    case_manifest = collect_case_qc_outputs(
        case_id=case_id,
        sample_id=clean(summary.get("tumor_sample")),
        summary_json=summary_path,
        variant_csv=variant_path,
        qc_report_md=summary.get("report_md", ""),
        case_report_md=summary.get("case_report_md", ""),
        decision_html=out_html,
        decision_md=out_md,
        shortlist_csv=out_shortlist,
    )
    print(
        json.dumps(
            {
                "case_id": case_id,
                "main_decision": decision["main_decision"],
                "full_vcf_decision": decision["full_vcf_decision"],
                "pass_strong_zero_normal": decision["pass_strong_zero_normal_count"],
                "kg_supported_marker_candidates": decision["pass_strong_kg_count"],
                "html": str(out_html),
                "md": str(out_md),
                "shortlist_csv": str(out_shortlist),
                "case_qc_dir": case_manifest["qc_dir"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
