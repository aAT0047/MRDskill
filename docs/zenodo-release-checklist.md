# Zenodo versioned archive checklist

No DOI is currently present in this checkout. A DOI must be copied from the actual Zenodo archive after a maintainer creates a release; do not guess, reserve, or fabricate one.

## Before release

- Confirm the canonical software title, author names, affiliations and ORCIDs with the maintainers. The GitHub login or license copyright line alone is not sufficient author metadata.
- Review the target version's source tree, tests, `LICENSE`, `docs/data_availability.md`, and all demo fixtures. Ensure no patient data, restricted article text, full indexes, secrets, local paths, or unlicensed benchmark data are included.
- Run `python tools/validate_skill_coverage.py` and `python scripts/run_demo.py --output-dir artifacts/demo --verify-expected` from a clean checkout.
- Review/update release notes, requirements, citation guidance, and the demo/data license statement.
- Add a `CITATION.cff` only after the canonical authors and release metadata are confirmed. Populate its DOI fields only after Zenodo mints the version DOI; distinguish the version-specific DOI from the concept DOI.

## Archive

1. A repository maintainer enables the GitHub integration in the intended Zenodo account and verifies the correct repository is selected.
2. Create a versioned GitHub Release from a reviewed, immutable tag (for example, the project's agreed `vX.Y.Z`). Do not reuse or move an already published tag.
3. Confirm Zenodo completed the archive and that its source snapshot, metadata, license, version, and downloadable files match the GitHub release.
4. Record both the version DOI and, if exposed, concept DOI from the Zenodo record. Add the version DOI to the release notes, README, and confirmed `CITATION.cff`; add the concept DOI only where it is clearly labeled.
5. Re-download or inspect the archived source snapshot and rerun the offline demo against that version. Save the verified archive URL and release/tag commit in project documentation.

Enabling Zenodo, creating a public GitHub release/tag, and minting a DOI are maintainer actions. This repository change prepares the instructions only; it does not publish or archive a release.
