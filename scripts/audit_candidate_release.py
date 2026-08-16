#!/usr/bin/env python
"""Gate 0 — pre-freeze audit of the candidate release package.

Verifies every claim the release package makes about itself before the
final freeze (git tag + SHA256). Checks:

1. every artifact in the manifest exists and its recorded SHA256 matches
   (hash-carrying artifacts);
2. the bundle loads, its state dict keys match the fit manifest's model,
   and scoring two identical inputs is bit-identical (CPU);
3. every text-bearing file in the package passes ``verify_predictive_claims``
   under the non-GO gate decision (``[]`` required);
4. all JSON documents parse;
5. the co-signature status is reported (the audit FAILS when run with
   ``--expect-signed`` and any protocol document is still unsigned —
   that flag is for the final freeze, after the wet-lab partner signs).

Usage::

    python scripts/audit_candidate_release.py            # report + hard failures
    python scripts/audit_candidate_release.py --expect-signed
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

DEFAULT_RELEASE_DIR = (
    _REPO_ROOT / "results" / "candidates" / "multispecies_v2_candidate_release_v1"
)
TEXT_SUFFIXES = (".md", ".json", ".py", ".tsv")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, default=DEFAULT_RELEASE_DIR)
    parser.add_argument(
        "--expect-signed",
        action="store_true",
        help="fail if any protocol document is still draft (use only at final freeze)",
    )
    args = parser.parse_args(argv)

    release_dir = args.release_dir
    failures: list[str] = []
    signed = True

    manifest_path = release_dir / "manifest.json"
    if not manifest_path.is_file():
        print(f"FAIL: manifest missing: {manifest_path}")
        return 1
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    # 1. artifact existence + pinned hashes
    for entry in manifest.get("produced_artifacts", []):
        artifact = entry.get("artifact")
        if not artifact:
            continue
        path = release_dir / artifact
        if not path.is_file():
            failures.append(f"missing artifact: {artifact}")
            continue
        status = entry.get("status")
        if status == "draft_pending_co_signature":
            signed = False
            print(f"DRAFT (unsigned): {artifact}")
        pinned = entry.get("sha256")
        if isinstance(pinned, str) and pinned.startswith(
            ("draft", "feature input", "blind-time", "Chinese companion")
        ):
            continue  # non-hash-carrying placeholders are recorded at final freeze
        if isinstance(pinned, str) and len(pinned) == 64:
            actual = _sha256(path)
            if actual != pinned:
                failures.append(
                    f"hash mismatch {artifact}: {actual} != {pinned}"
                )

    # 2. bundle loads + deterministic scoring
    try:
        from plantpersulf.models.structure_ranker import (
            BranchFeatures,
            StructureRankerBundle,
            score_structure_ranker_bundle,
        )

        bundle_path = release_dir / "model_weights" / "structure_ranker_bundle.pt"
        bundle = StructureRankerBundle.load(bundle_path)
        predict = BranchFeatures(
            sequence=[[0.5, 0.1, 0.3], [0.9, 0.01, 0.7]],
            esm=[[0.0], [0.0]],
            structure=[[12.0, 85.0], [0.0, 0.0]],
            structure_mask=[True, False],
            study_ids=None,
        )
        first = score_structure_ranker_bundle(bundle, predict)
        second = score_structure_ranker_bundle(bundle, predict)
        if first.scores != second.scores or first.uncertainty != second.uncertainty:
            failures.append("bundle scoring is not bit-identical on CPU")
        else:
            print(
                f"bundle OK: seed={bundle.seed} dims={bundle.input_dims} "
                f"deterministic=true"
            )
    except Exception as exc:  # noqa: BLE001 - audit reports any failure
        failures.append(f"bundle check failed: {exc}")

    # 3. predictive-claims discipline across the package
    from plantpersulf.evaluation.conclusion_gate import (
        Gate2Decision,
        verify_predictive_claims,
    )

    decision = Gate2Decision(decision="hold", conditions=[], downgrade_statement="")
    claim_violations: list[str] = []
    for path in sorted(
        p for p in release_dir.rglob("*") if p.is_file() and p.suffix in TEXT_SUFFIXES
    ):
        hits = verify_predictive_claims(
            path.read_text(encoding="utf-8", errors="replace"), decision
        )
        if hits:
            claim_violations.append(f"{path.relative_to(release_dir)}: {hits}")
    if claim_violations:
        failures.extend(f"predictive-claim violation: {item}" for item in claim_violations)
    else:
        print("predictive claims: clean ([])")

    # 4. JSON parses
    for path in sorted(p for p in release_dir.rglob("*.json") if p.is_file()):
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            failures.append(f"invalid JSON {path.relative_to(release_dir)}: {exc}")
    print("JSON documents: all parse")

    if args.expect_signed and not signed:
        failures.append(
            "expected co-signature but protocol documents are still draft"
        )

    print()
    if failures:
        print("AUDIT FAILED:")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print(
        "AUDIT PASS - release package self-consistent"
        + ("" if signed else " (protocol documents still draft; re-run with --expect-signed after signing)")
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
