"""The software is owned by the author alone; citations carry the SPL/VSSC/ISRO
affiliation but not the visitor title."""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
# The installed-package CI job runs the tests outside the source tree.
pytestmark = pytest.mark.skipif(not (ROOT / "LICENSE").exists(), reason="source tree not present")


def test_license_names_no_institution():
    lic = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert "Copyright (c) 2026 Keshav Aggarwal" in lic
    for banned in ("ISRO", "VSSC", "Space Physics Laboratory"):
        assert banned not in lic, banned


def test_citation_affiliation_without_visitor_title():
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert 'affiliation: "Space Physics Laboratory (SPL), Vikram Sarabhai Space Centre (VSSC), ISRO"' in cff
    assert "Student Visitor" not in cff
