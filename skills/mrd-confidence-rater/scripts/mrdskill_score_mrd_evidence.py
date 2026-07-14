#!/usr/bin/env python3
"""Rule-score MRD evidence candidates for graph readiness."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


ADDED_COLUMNS = [
    "evidence_source_type_final",
    "confidence_score",
    "confidence_reason",
    "evidence_strength_score",
    "graph_inclusion_status",
    "audit_priority",
    "downgrade_reason",
    "upgrade_reason",
    "contradiction_status",
    "reproducibility_status",
    "clinical_relevance_status",
    "biological_plausibility_status",
    "source_traceability_status",
    "final_review_note",
]

HIGH_SOURCE_TYPES = {"clinical_validation", "primary_result", "supplementary_result", "meta_analysis"}
WEAK_SOURCE_TYPES = {
    "abstract_only_claim",
    "discussion_or_background_claim",
    "review_statement",
    "inferred_candidate",
    "uncertain",
}

DIRECT_RELEVANCE = {
    "recurrence_marker",
    "relapse_driver",
    "prognostic_marker",
    "treatment_monitoring_marker",
    "tumor_specific_marker",
    "methylation_marker",
    "fragmentomics_marker",
}

CHIP_GENES = {
    "DNMT3A",
    "TET2",
    "ASXL1",
    "JAK2",
    "TP53",
    "PPM1D",
    "SF3B1",
    "SRSF2",
    "U2AF1",
    "IDH1",
    "IDH2",
    "ATM",
    "CHEK2",
    "GNAS",
    "KRAS",
}

VARIANT_CONTEXT_RE = re.compile(
    r"\bvariant\b|\bmutation\b|\bmutant\b|\bVAF\b|variant allele|allele frequency|"
    r"\bAF\b|\bctDNA\b|\bcfDNA\b|\bplasma\b|\btumou?r\b|\btissue\b|\bsample\b|"
    r"baseline|postoperative|post[- ]?treatment|timepoint|MRD",
    re.IGNORECASE,
)


def is_reported(value: str | None) -> bool:
    return bool(value and value not in {"not_reported", "unavailable", "uncertain", "not_assessable", ""})


def split_genes(value: str) -> set[str]:
    return {item.strip().upper() for item in (value or "").split(";") if item.strip() and item != "not_reported"}


def is_supplementary(row: dict[str, str]) -> bool:
    return any(
        is_reported(row.get(key))
        for key in ["supplementary_file", "supplementary_table", "source_sheet", "source_row"]
    ) or row.get("source_section") in {"Supplementary_Table", "Supplementary_Figure", "Supplementary_Dataset"}


def has_structured_variant_evidence(row: dict[str, str]) -> bool:
    text = " ".join(
        [
            row.get("gene", ""),
            row.get("variant", ""),
            row.get("variant_type", ""),
            row.get("quoted_or_paraphrased_evidence", ""),
        ]
    )
    has_gene = is_reported(row.get("gene"))
    has_variant = is_reported(row.get("variant")) or row.get("variant_type") in {
        "SNV",
        "indel",
        "fusion",
        "CNV",
        "SV",
        "protein_change",
        "methylation",
        "fragmentomics",
    }
    return (has_gene and has_variant) or (has_gene and bool(VARIANT_CONTEXT_RE.search(text)))


def has_statistical_or_validation_support(row: dict[str, str], upgrades: list[str]) -> bool:
    return any(
        reason in upgrades
        for reason in [
            "strong_statistics",
            "validation_cohort",
            "prospective_study",
            "survival_endpoint",
            "clinically_validated_assay",
        ]
    ) or is_reported(row.get("sample_size"))


def source_traceability(row: dict[str, str]) -> str:
    has_paper = is_reported(row.get("doi")) or is_reported(row.get("pmid")) or is_reported(row.get("pmcid"))
    has_section = is_reported(row.get("source_section")) and row.get("source_section") != "unknown"
    has_quote = is_reported(row.get("quoted_or_paraphrased_evidence"))
    has_location = any(
        is_reported(row.get(key))
        for key in ["source_table", "source_figure", "supplementary_file", "supplementary_table", "source_sheet", "source_row"]
    )
    if has_paper and has_section and (has_quote or has_location):
        return "fully_traceable"
    if has_paper and (has_section or has_quote):
        return "partially_traceable"
    if has_paper:
        return "weakly_traceable"
    return "not_traceable"


def clinical_relevance(row: dict[str, str]) -> str:
    if row.get("mrd_relevance") in DIRECT_RELEVANCE:
        return "directly_MRD_related"
    if any(is_reported(row.get(key)) for key in ["recurrence_association", "relapse_association", "prognosis_association"]):
        return "directly_MRD_related"
    if is_reported(row.get("assay_type")) or "ctDNA" in row.get("quoted_or_paraphrased_evidence", ""):
        return "indirectly_MRD_related"
    if is_supplementary(row) and has_structured_variant_evidence(row):
        return "structured_variant_evidence"
    return "not_directly_MRD_related"


def biological_plausibility(row: dict[str, str]) -> str:
    genes = split_genes(row.get("gene_symbol_normalized", ""))
    text = " ".join(
        [
            row.get("biomarker_type", ""),
            row.get("mrd_relevance", ""),
            row.get("quoted_or_paraphrased_evidence", ""),
        ]
    ).lower()
    if "chip" in text or genes & CHIP_GENES and row.get("matched_WBC_used") != "yes":
        return "likely_CHIP" if "matched wbc" in text or "clonal hematopoiesis" in text else "biologically_uncertain"
    if row.get("mrd_relevance") in DIRECT_RELEVANCE and is_reported(row.get("cancer_type")):
        return "biologically_plausible"
    if is_reported(row.get("gene")) or is_reported(row.get("assay_type")):
        return "biologically_uncertain"
    return "not_assessable"


def upgrade_reasons(row: dict[str, str]) -> list[str]:
    reasons: list[str] = []
    text = " ".join([row.get("cohort_type", ""), row.get("validation_status", ""), row.get("quoted_or_paraphrased_evidence", "")]).lower()
    if "prospective" in text:
        reasons.append("prospective_study")
    if "validation" in text or row.get("validation_status") in {"validated", "clinical_study"}:
        reasons.append("validation_cohort")
    if is_reported(row.get("hazard_ratio")) or is_reported(row.get("p_value")):
        reasons.append("strong_statistics")
    if any(is_reported(row.get(key)) for key in ["source_table", "source_figure"]):
        reasons.append("main_table_or_figure_support")
    if is_supplementary(row):
        reasons.append("supplementary_table_support")
    if has_structured_variant_evidence(row):
        reasons.append("structured_variant_evidence")
    if row.get("source_section") == "Results":
        reasons.append("main_table_or_figure_support")
    if row.get("matched_WBC_used") == "yes":
        reasons.append("matched_WBC_filtering")
    if row.get("tumor_informed_or_tumor_naive") == "tumor_informed":
        reasons.append("tumor_tissue_confirmation")
    if any(row.get(key) == "reported" for key in ["recurrence_association", "relapse_association"]):
        reasons.append("direct_recurrence_endpoint")
    if is_reported(row.get("survival_endpoint")):
        reasons.append("survival_endpoint")
    return sorted(set(reasons))


def downgrade_reasons(row: dict[str, str], traceability: str) -> list[str]:
    reasons: list[str] = []
    source_type = row.get("evidence_source_type_preliminary", "uncertain")
    if source_type == "abstract_only_claim":
        reasons.append("abstract_only")
    if source_type == "discussion_or_background_claim":
        reasons.append("discussion_only")
    if source_type == "review_statement":
        reasons.append("review_only")
    structured_variant = has_structured_variant_evidence(row)
    if not structured_variant and not is_reported(row.get("hazard_ratio")) and not is_reported(row.get("p_value")):
        reasons.append("missing_statistics")
    if row.get("validation_status") in {"not_reported", "", None}:
        reasons.append("no_validation_cohort")
    if not is_reported(row.get("cancer_type")):
        reasons.append("unclear_cancer_context")
    if not is_reported(row.get("gene")) and row.get("variant_type") not in {"methylation", "fragmentomics"}:
        reasons.append("unclear_gene_context")
    if not is_reported(row.get("assay_type")):
        reasons.append("unclear_assay_context")
    if row.get("matched_WBC_used") != "yes":
        reasons.append("no_matched_WBC")
    genes = split_genes(row.get("gene_symbol_normalized", ""))
    if genes & CHIP_GENES and row.get("matched_WBC_used") != "yes":
        reasons.append("possible_CHIP")
    if traceability in {"weakly_traceable", "not_traceable"}:
        reasons.append("weak_traceability")
    if clinical_relevance(row) == "not_directly_MRD_related":
        reasons.append("not_directly_MRD_related")
    return sorted(set(reasons))


def score_row(row: dict[str, str]) -> dict[str, str]:
    final_type = row.get("evidence_source_type_preliminary") or "uncertain"
    traceability = source_traceability(row)
    relevance = clinical_relevance(row)
    plausibility = biological_plausibility(row)
    upgrades = upgrade_reasons(row)
    downgrades = downgrade_reasons(row, traceability)

    strength = 1
    if final_type in HIGH_SOURCE_TYPES:
        strength += 1
    if relevance == "directly_MRD_related":
        strength += 1
    if relevance == "structured_variant_evidence":
        strength += 1
    if traceability in {"fully_traceable", "partially_traceable"}:
        strength += 1
    if "strong_statistics" in upgrades or "validation_cohort" in upgrades:
        strength += 1
    if "supplementary_table_support" in upgrades and "structured_variant_evidence" in upgrades:
        strength += 1
    if final_type in WEAK_SOURCE_TYPES:
        strength -= 1
    if "possible_CHIP" in downgrades or "weak_traceability" in downgrades:
        strength -= 1
    strength = max(0, min(5, strength))

    has_support = has_statistical_or_validation_support(row, upgrades)

    if traceability == "not_traceable" or relevance == "not_directly_MRD_related":
        confidence = "Exclude" if strength == 0 else "Low"
    elif strength >= 4 and final_type in HIGH_SOURCE_TYPES and relevance == "directly_MRD_related" and has_support:
        confidence = "High"
    elif strength >= 4 and final_type == "supplementary_result" and relevance == "structured_variant_evidence":
        confidence = "Medium"
    elif strength >= 4 and final_type == "supplementary_result" and not has_support:
        confidence = "Medium"
    elif strength >= 3:
        confidence = "Medium"
    elif strength >= 1:
        confidence = "Low"
    else:
        confidence = "Exclude"

    if confidence == "High" and strength >= 4:
        inclusion = "auto_include_high_confidence"
    elif relevance == "structured_variant_evidence" and confidence in {"Medium", "Low"}:
        inclusion = "include_as_medium_confidence" if confidence == "Medium" else "candidate_only"
    elif confidence == "Medium":
        inclusion = "include_as_medium_confidence"
    elif confidence == "Low":
        inclusion = "include_as_low_confidence_annotation"
    elif confidence == "Exclude":
        inclusion = "exclude"
    else:
        inclusion = "candidate_only"

    audit_priority = "no_review_needed"
    if "possible_CHIP" in downgrades or plausibility in {"likely_CHIP", "biologically_uncertain"}:
        audit_priority = "urgent_review"
    elif relevance == "structured_variant_evidence":
        audit_priority = "high_review"
    elif is_supplementary(row) and confidence == "Medium":
        audit_priority = "medium_review"
    elif confidence in {"Low", "Uncertain"} or traceability != "fully_traceable":
        audit_priority = "high_review"
    elif confidence == "Medium":
        audit_priority = "medium_review"

    confidence_reason = (
        f"Rule score {strength}/5 from source_type={final_type}, "
        f"traceability={traceability}, clinical_relevance={relevance}."
    )

    row.update(
        {
            "evidence_source_type_final": final_type,
            "confidence_score": confidence,
            "confidence_reason": confidence_reason,
            "evidence_strength_score": str(strength),
            "graph_inclusion_status": inclusion,
            "audit_priority": audit_priority,
            "downgrade_reason": ";".join(downgrades) if downgrades else "not_reported",
            "upgrade_reason": ";".join(upgrades) if upgrades else "not_reported",
            "contradiction_status": "not_assessable",
            "reproducibility_status": "not_assessable",
            "clinical_relevance_status": relevance,
            "biological_plausibility_status": plausibility,
            "source_traceability_status": traceability,
            "final_review_note": "python_rule_scored; review before clinical use",
        }
    )
    return row


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    return columns, rows


def write_csv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_jsonl(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Rule-score raw MRD evidence candidates.")
    parser.add_argument("--input-csv", default="data/evidence/mrd_evidence_raw.csv")
    parser.add_argument("--out-csv", default="data/evidence/mrd_evidence_scored.csv")
    parser.add_argument("--out-jsonl", default="data/evidence/mrd_evidence_scored.jsonl")
    args = parser.parse_args()

    columns, rows = read_csv(Path(args.input_csv))
    scored = [score_row(row) for row in rows]
    out_columns = columns + [column for column in ADDED_COLUMNS if column not in columns]
    write_csv(Path(args.out_csv), out_columns, scored)
    write_jsonl(Path(args.out_jsonl), scored)
    print(f"Wrote scored evidence CSV: {Path(args.out_csv).resolve()} ({len(scored)} rows)")
    print(f"Wrote scored evidence JSONL: {Path(args.out_jsonl).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
