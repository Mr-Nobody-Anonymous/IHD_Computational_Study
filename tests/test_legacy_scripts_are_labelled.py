"""The legacy root-level scripts must stay clearly marked as non-usable.

They were the target of an external critique that assumed they were the study. They
are not, and nothing in the package imports them - but a reader who finds one has to be
told that before reading a number. These tests fail the build if that labelling is
removed or if a legacy script is wired into the tested package.
"""

from __future__ import annotations

import ast
import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
SRC = REPO / "src"
LEGACY = sorted(p for p in SRC.glob("0*.py") if p.name[0].isdigit())
PACKAGE = SRC / "ihd_solutions"


def test_the_five_legacy_scripts_are_still_present():
    assert [p.name for p in LEGACY] == [
        "01_multistructure_sasa_benchmark.py",
        "02_interface_contact_mapping.py",
        "03_plaque_biomechanics_model.py",
        "04_immune_cytokine_network.py",
        "05_crispr_genomic_safety_cfd.py",
    ]


@pytest.mark.parametrize("path", LEGACY, ids=lambda p: p.name)
def test_each_legacy_script_declares_itself_non_usable(path):
    docstring = ast.get_docstring(ast.parse(path.read_text())) or ""
    assert "NON-USABLE" in docstring
    assert "What it can support:" in docstring
    assert "What it CANNOT support:" in docstring


@pytest.mark.parametrize("path", LEGACY, ids=lambda p: p.name)
def test_each_legacy_script_warns_when_run(path):
    source = path.read_text()
    assert "_warnings.warn(" in source
    assert "NOT findings" in source


@pytest.mark.parametrize("path", LEGACY, ids=lambda p: p.name)
def test_each_legacy_script_parses(path):
    ast.parse(path.read_text())


def test_the_plaque_script_admits_its_invented_mmp_activities():
    """The four regimen multipliers have no source. The banner must say so."""
    text = (SRC / "03_plaque_biomechanics_model.py").read_text()
    assert "NO source" in text
    assert "HYPOTHETICAL" in text


def test_the_crispr_script_does_not_claim_a_genome_wide_screen():
    text = (SRC / "05_crispr_genomic_safety_cfd.py").read_text()
    assert "CANDIDATE-LIST" in text
    assert "not a genome" in text
    docstring = ast.get_docstring(ast.parse(text)) or ""
    # The file name says "genomic_safety"; the docstring may not promise genome-wide.
    assert "genome-wide off-target safety assessment" in docstring


def test_no_package_module_imports_a_legacy_script():
    """If the package ever starts depending on them, the label becomes a lie."""
    pattern = re.compile(r"^\s*(?:from|import)\s+0\d_", re.MULTILINE)
    for module in PACKAGE.glob("*.py"):
        assert not pattern.search(module.read_text()), module.name


def test_the_package_itself_never_silences_the_warning():
    for module in PACKAGE.glob("*.py"):
        assert "filterwarnings" not in module.read_text(), module.name
