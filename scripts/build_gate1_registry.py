#!/usr/bin/env python
"""Build the gate-1 (oxidizable-Cys) site registry TSV.

Stage-1 training labels for the two-step persulfidation model: cysteine
sites known to be oxidized (S-nitrosylation first — the regulatory-oxidation
class best aligned with the review's "oxiPTMs compete for the same reactive
Cys" framing). Coordinate verification is fail-closed against the registered
reference proteomes; a site is ``verified`` only when its accession is in a
registered proteome AND the reported position is a Cys.

Usage::

    python scripts/build_gate1_registry.py \\
        --output data/registry/gate1_oxidation_sites_v1.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from plantpersulf.evidence.gate1_oxidation import TIER1_SNO, gate1_audit  # noqa: E402
from plantpersulf.evidence.oxiptm_sites import (  # noqa: E402
    ARABIDOPSIS,
    RICE,
    TOMATO,
)
from plantpersulf.features.sequence import _load_proteome  # noqa: E402

OUTPUT = _REPO_ROOT / "data" / "registry" / "gate1_oxidation_sites_v1.tsv"
# Tier-2: 2935 verified Arabidopsis SNO sites extracted from the OA workflow
# paper PMC12454217 (MOESM2), TAIR loci mapped to UniProt and coordinate-verified.
ARABI_SNO_VERIFIED = (
    _REPO_ROOT
    / "data"
    / "raw"
    / "supplements"
    / "PXD_tomato_SNO"
    / "arabi_sno_sites_verified.tsv"
)

PROTEOME_PATHS = {
    ARABIDOPSIS: (
        _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
    ),
    TOMATO: (
        _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
    ),
    RICE: (
        _REPO_ROOT / "data" / "raw" / "references" / "rice_proteome_v1"
        / "uniprot_rice_v1.fasta"
    ),
}


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    return parser


def _load_arabi_sno_verified(path: Path) -> list[dict[str, object]]:
    """Tier-2 bulk SNO sites already mapped + coordinate-verified."""
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(
                {
                    "gene": row["tair_locus"],
                    "species": ARABIDOPSIS,
                    "cys_position": int(row["cys_position"]),
                    "oxidation_type": "s-nitrosylation",
                    "accession": row["uniprot_accession"],
                    "observed_residue": "C",
                    "status": "verified",
                    "source": "PMC12454217 MOESM2 (S-nitrosylation column)",
                    "source_note": (
                        "Novel proteomics workflow paper; TAIR locus -> UniProt "
                        "mapped, coordinate-verified in registered proteome v2"
                    ),
                }
            )
    return rows


def main(argv: list[str] | None = None) -> None:
    args = build_argument_parser().parse_args(argv)
    proteomes = {
        species: _load_proteome(path) for species, path in PROTEOME_PATHS.items()
    }
    rows = gate1_audit(TIER1_SNO, proteomes=proteomes)
    rows.extend(_load_arabi_sno_verified(ARABI_SNO_VERIFIED))
    fieldnames = [
        "gene",
        "species",
        "cys_position",
        "oxidation_type",
        "accession",
        "observed_residue",
        "status",
        "source",
        "source_note",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    n_verified = sum(1 for r in rows if r["status"] == "verified")
    print(
        f"gate-1 registry: {len(rows)} sites "
        f"({n_verified} verified in registered proteomes)"
    )
    for r in rows:
        print(
            f"  {r['status']:<13} {r['gene']:<8} {r['species']:<11} "
            f"C{r['cys_position']} {r['accession'] or '-'}"
        )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
