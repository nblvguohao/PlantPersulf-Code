"""Fail-closed scientific input audits."""

import csv
import hashlib
from pathlib import Path
from typing import Any

import yaml


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_registered_input(input_path: Path, registry_path: Path) -> None:
    """Reject an input unless its path and current SHA256 are registered."""
    input_path = input_path.resolve()
    registry_path = registry_path.resolve()
    if not registry_path.is_file():
        raise RuntimeError(
            f"input is not registered: registry missing: {registry_path}"
        )

    with registry_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"path", "sha256"}.issubset(
            reader.fieldnames
        ):
            raise RuntimeError(
                "input is not registered: registry requires path and sha256"
            )
        for row in reader:
            registered_path = Path(row["path"])
            if not registered_path.is_absolute():
                registered_path = registry_path.parent / registered_path
            if registered_path.resolve() != input_path:
                continue
            actual_sha256 = _sha256(input_path)
            if actual_sha256.lower() != row["sha256"].strip().lower():
                raise RuntimeError(
                    f"SHA256 mismatch for registered input: {input_path}"
                )
            return

    raise RuntimeError(f"input is not registered: {input_path}")


def _forbidden_filename_terms(config_path: Path) -> tuple[str, ...]:
    if not config_path.is_file():
        raise RuntimeError(f"scientific integrity config missing: {config_path}")
    loaded: Any = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict):
        raise RuntimeError("scientific integrity config must be a mapping")
    release = loaded.get("release")
    if not isinstance(release, dict):
        raise RuntimeError("scientific integrity config requires release policy")
    terms = release.get("forbidden_filename_terms")
    if not isinstance(terms, list) or not terms:
        raise RuntimeError("release policy requires forbidden_filename_terms")
    if any(not isinstance(term, str) or not term.strip() for term in terms):
        raise RuntimeError("forbidden_filename_terms must contain non-empty strings")
    return tuple(term.strip().lower() for term in terms)


def assert_no_forbidden_scientific_artifacts(
    release_path: Path,
    config_path: Path,
) -> None:
    """Reject release files whose paths contain policy-forbidden terms."""
    if not release_path.is_dir():
        raise RuntimeError(f"release directory missing: {release_path}")
    forbidden_terms = _forbidden_filename_terms(config_path)
    for artifact in release_path.rglob("*"):
        if not artifact.is_file():
            continue
        relative_name = artifact.relative_to(release_path).as_posix().lower()
        matched_term = next(
            (term for term in forbidden_terms if term in relative_name),
            None,
        )
        if matched_term is not None:
            raise RuntimeError(
                "forbidden scientific artifact "
                f"'{artifact}' matches policy term '{matched_term}'"
            )
