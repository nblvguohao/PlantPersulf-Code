"""Traceable provenance records for scientific datasets."""

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class DatasetRecord:
    """Minimum provenance required before a dataset may be used."""

    accession: str
    repository: str
    source_url: str
    scientific_role: str
    metadata_retrieved_at: str
    metadata_sha256: str

    def __post_init__(self) -> None:
        required_fields = (
            "accession",
            "repository",
            "source_url",
            "scientific_role",
            "metadata_retrieved_at",
            "metadata_sha256",
        )
        for field_name in required_fields:
            if not getattr(self, field_name).strip():
                raise ValueError(f"{field_name} is required")
        if re.fullmatch(r"[0-9a-fA-F]{64}", self.metadata_sha256) is None:
            raise ValueError(
                "metadata_sha256 must be a 64-character hexadecimal digest"
            )
