#!/usr/bin/env python
"""Build the wet-lab handoff package for the 600-site blind cohort.

Reads the FROZEN release artifacts (top_k_candidates.tsv,
matched_controls.tsv) and the SHA-pinned tomato reference proteome, then
writes a blinded assay list for the Zhang Hua laboratory:

- ``assay_list_blind.tsv`` — 600 sites (Top-200 candidates + 2 nearest
  matched controls each), one row per site, identified ONLY by blind_id
  (BLIND-0001..BLIND-0600), deterministically shuffled (seed 20260813).
  This is the ONLY site-level file sent to the wet lab.
- ``blind_id_key.DO_NOT_SEND.tsv`` — blind_id -> role/score mapping.
  Stays with the modeling side; its SHA256 is recorded in the handoff
  manifest so the key can be verified at unblinding.
- ``handoff_manifest.json`` — input/output hashes, shuffle seed, counts.

Role-blinding rationale: the wet laboratory must not know which sites are
candidates and which are matched controls, otherwise assay QC could
differ systematically between the two groups. The modeling side must not
see labels until the pre-registered unblinding step (SAP protocol).

This script changes no frozen artifact; it only derives the handoff view.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

RELEASE_DIR = _REPO_ROOT / "results" / "candidates" / "multispecies_v2_candidate_release_v1"
PROTEOME = _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
DEFAULT_OUTPUT_DIR = _REPO_ROOT / "results" / "handoff" / "zhang_lab_blind_cohort_v1"

K_CANDIDATES = 200
CONTROLS_PER_CANDIDATE = 2
WINDOW = 10  # +/- residues around the target Cys
SHUFFLE_SEED = 20260813


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_proteome(path: Path) -> dict[str, dict[str, str]]:
    """accession -> {header, gene, protein_name, sequence}."""
    entries: dict[str, dict[str, str]] = {}
    header = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                _store(entries, header, chunks)
            header = line[1:]
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        _store(entries, header, chunks)
    return entries


def _store(entries: dict[str, dict[str, str]], header: str, chunks: list[str]) -> None:
    parts = header.split("|")
    accession = parts[1] if len(parts) >= 3 else header.split()[0]
    gene_match = re.search(r"GN=(\S+)", header)
    name_match = re.search(r"^\S+\s+(.*?)\s+OS=", header)
    entries[accession] = {
        "header": header,
        "gene": gene_match.group(1) if gene_match else "",
        "protein_name": name_match.group(1) if name_match else "",
        "sequence": "".join(chunks),
    }


def _window(sequence: str, position: int) -> str:
    """±WINDOW residues around 1-based Cys ``position``; target shown as [C]."""
    idx = position - 1
    lo = max(0, idx - WINDOW)
    hi = min(len(sequence), idx + WINDOW + 1)
    return sequence[lo:idx] + "[C]" + sequence[idx + 1:hi]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(argv)

    candidates_path = RELEASE_DIR / "top_k_candidates.tsv"
    controls_path = RELEASE_DIR / "matched_controls.tsv"
    for path in (candidates_path, controls_path, PROTEOME):
        if not path.is_file():
            raise RuntimeError(f"missing input: {path}")

    with candidates_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    head = rows[:K_CANDIDATES]
    if len(head) < K_CANDIDATES:
        raise RuntimeError(f"candidate table shorter than K={K_CANDIDATES}")

    controls_by_candidate: dict[tuple[str, str], list[dict[str, str]]] = {}
    with controls_path.open(encoding="utf-8", newline="") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            key = (record["candidate_site_key"], record["candidate_cys_position"])
            controls_by_candidate.setdefault(key, []).append(record)

    proteome = _parse_proteome(PROTEOME)
    print(f"proteome entries: {len(proteome)}")

    sites: list[dict[str, str]] = []  # pre-shuffle rows with role info
    shortfalls: list[str] = []
    for rank, cand in enumerate(head, start=1):
        key = (cand["site_key"], cand["cys_position"])
        controls = controls_by_candidate.get(key, [])
        # matched_controls.tsv is written nearest-first; take the 2 nearest
        chosen = controls[:CONTROLS_PER_CANDIDATE]
        if len(chosen) < CONTROLS_PER_CANDIDATE:
            shortfalls.append(
                f"{cand['site_key']}:{cand['cys_position']} has only {len(chosen)} controls"
            )
        sites.append(
            {
                "role": "candidate",
                "rank": str(rank),
                "site_key": cand["site_key"],
                "cys_position": cand["cys_position"],
                "score": cand["score"],
                "uncertainty": cand["uncertainty"],
                "pair_rank": "",
            }
        )
        for control_rank, control in enumerate(chosen, start=1):
            sites.append(
                {
                    "role": f"control_{control_rank}",
                    "rank": str(rank),
                    "site_key": control["control_site_key"],
                    "cys_position": control["control_cys_position"],
                    "score": "",
                    "uncertainty": "",
                    "pair_rank": str(control_rank),
                }
            )
    n_candidates = sum(1 for s in sites if s["role"] == "candidate")
    n_controls = len(sites) - n_candidates
    print(f"cohort: {n_candidates} candidates + {n_controls} controls = {len(sites)} sites")
    if shortfalls:
        print(f"WARNING: {len(shortfalls)} candidates with fewer than 2 controls:")
        for item in shortfalls:
            print(f"  {item}")

    # --- annotate from proteome ------------------------------------------------
    missing: list[str] = []
    for site in sites:
        accession = site["site_key"].split("|", 1)[1]
        entry = proteome.get(accession)
        position = int(site["cys_position"])
        if entry is None or position > len(entry["sequence"]):
            missing.append(site["site_key"])
            continue
        if entry["sequence"][position - 1] != "C":
            raise RuntimeError(f"coordinate check failed: {site['site_key']}:{position}")
        site["accession"] = accession
        site["gene"] = entry["gene"]
        site["protein_name"] = entry["protein_name"]
        site["protein_length"] = str(len(entry["sequence"]))
        site["sequence_window"] = _window(entry["sequence"], position)
    if missing:
        raise RuntimeError(f"{len(missing)} sites missing from proteome: {missing[:5]}")

    # --- deterministic shuffle + blind ids -------------------------------------
    order = list(range(len(sites)))
    random.Random(SHUFFLE_SEED).shuffle(order)
    for blind_index, site_index in enumerate(order, start=1):
        sites[site_index]["blind_id"] = f"BLIND-{blind_index:04d}"

    args.output_dir.mkdir(parents=True, exist_ok=True)

    blind_path = args.output_dir / "assay_list_blind.tsv"
    with blind_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "blind_id", "accession", "gene", "protein_name",
                "cys_position", "protein_length", "sequence_window",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for site_index in order:
            site = sites[site_index]
            writer.writerow(
                {
                    "blind_id": site["blind_id"],
                    "accession": site["accession"],
                    "gene": site["gene"],
                    "protein_name": site["protein_name"],
                    "cys_position": site["cys_position"],
                    "protein_length": site["protein_length"],
                    "sequence_window": site["sequence_window"],
                }
            )

    key_path = args.output_dir / "blind_id_key.DO_NOT_SEND.tsv"
    with key_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "blind_id", "role", "candidate_rank", "site_key", "cys_position",
                "score", "uncertainty", "control_pair_rank",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for site_index in order:
            site = sites[site_index]
            writer.writerow(
                {
                    "blind_id": site["blind_id"],
                    "role": site["role"],
                    "candidate_rank": site["rank"],
                    "site_key": site["site_key"],
                    "cys_position": site["cys_position"],
                    "score": site["score"],
                    "uncertainty": site["uncertainty"],
                    "control_pair_rank": site["pair_rank"],
                }
            )

    manifest = {
        "handoff_id": "zhang_lab_blind_cohort_v1",
        "release_id": "multispecies-v2-candidate-release-v1",
        "design": {
            "k_candidates": K_CANDIDATES,
            "controls_per_candidate": CONTROLS_PER_CANDIDATE,
            "total_sites": len(sites),
            "n_candidates": n_candidates,
            "n_controls": n_controls,
            "control_shortfalls": shortfalls,
        },
        "blinding": {
            "shuffle_seed": SHUFFLE_SEED,
            "rule": "wet laboratory receives assay_list_blind.tsv only (blind_id, no role, no score); "
                    "blind_id_key.DO_NOT_SEND.tsv stays with the modeling side and its SHA256 is "
                    "recorded here; the key is verified against this hash at unblinding",
        },
        "inputs": {
            "candidates": {"path": str(candidates_path), "sha256": _sha256(candidates_path)},
            "matched_controls": {"path": str(controls_path), "sha256": _sha256(controls_path)},
            "proteome": {"path": str(PROTEOME), "sha256": _sha256(PROTEOME)},
        },
        "outputs": {
            "assay_list_blind": {"file": blind_path.name, "sha256": _sha256(blind_path)},
            "blind_id_key": {"file": key_path.name, "sha256": _sha256(key_path)},
        },
        "created": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path = args.output_dir / "handoff_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"wrote {blind_path}")
    print(f"wrote {key_path}")
    print(f"wrote {manifest_path}")


if __name__ == "__main__":
    main()
