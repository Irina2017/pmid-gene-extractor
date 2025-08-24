"""
Test CLI functionality
"""

import time

import pytest
from click.testing import CliRunner

from pmid_parser.config import RATE_LIMIT_DELAY_SECONDS
from pmid_parser.main import main


@pytest.fixture
def cli_runner():
    """Provide a Click CLI test runner."""
    return CliRunner()


# to avoid API problems with PubMed or PMC
@pytest.fixture(autouse=True)
def rate_guard():
    time.sleep(RATE_LIMIT_DELAY_SECONDS)


def test_cli_requires_pmid(cli_runner):
    """Test that CLI requires PMID argument."""
    result = cli_runner.invoke(main, [])
    assert result.exit_code != 0
    assert "Missing option '--pmid'" in result.output


def test_cli_works(cli_runner):
    """Test CLI accepts PMID and doesn't crash."""
    result = cli_runner.invoke(main, ["--pmid", "38790019"])
    assert result.exit_code == 0


@pytest.mark.parametrize(
    "pmid_input",
    [
        "38790019",
        "  38790019  ",  # With whitespace
        "PMID: 38790019",  # With prefix
        "  PMID:    38790019  ",  # With prefix and whitespace
    ],
)
def test_cli_accepts_valid_pmid(cli_runner, pmid_input):
    """Test CLI accepts and cleans valid PMID formats."""

    result = cli_runner.invoke(main, ["--pmid", pmid_input])
    assert result.exit_code == 0


def test_cli_run_test(cli_runner):
    """Test CLI accepts and run test PMID - integrated test"""

    result = cli_runner.invoke(main, ["--pmid", "38790019"])
    assert result.exit_code == 0
    assert result.output.splitlines()[0].strip() == "Processing PMID: 38790019"
    expected = [
        # --- ANALYZING lines ---
        "ANALYZING APOL1 - Original mentions: 3",
        "ANALYZING RRAGD - Original mentions: 5",
        "ANALYZING HNF1A - Original mentions: 4",
        "ANALYZING COL4A3 - Original mentions: 6",
        "ANALYZING NPHS2 - Original mentions: 7",
        "ANALYZING HERC2 - Original mentions: 1",
        # --- Validated gene results header ---
        "Found 6 validated genes",
        # --- Check a few representative validated genes ---
        #    "1. APOL1 (HGNC:618)"
        "   Symbol: APOL1",
        "   Full name: apolipoprotein L1",
        "   Aliases: ",
        "   HG38 coordinates: 22:36253071-36267530",
        "   HG19 coordinates: 22:36649056-36663576",
        # "2. RRAGD (HGNC:19903)"
        "   Symbol: RRAGD",
        "   Full name: Ras related GTP binding D",
        "   Aliases: DKFZP761H171, bA11D8.2.1",
        "   HG38 coordinates: 6:89364616-89412273",
        "   HG19 coordinates: 6:90074355-90121989",
        # "3. HNF1A (HGNC:11621)"
        "   Symbol: HNF1A",
        "   Full name: HNF1 homeobox A",
        "   Aliases: HNF1, LFB1, HNF1α",
        "   HG38 coordinates: 12:120978543-121002512",
        "   HG19 coordinates: 12:121416346-121440315",
        # "4. COL4A3 (HGNC:2204)"
        "   Symbol: COL4A3",
        "   Full name: collagen type IV alpha 3 chain",
        "   Aliases: ",
        "   HG38 coordinates: 2:227164624-227314792",
        "   HG19 coordinates: 2:228029281-228179508",
        # "5. NPHS2 (HGNC:13394)"
        "   Symbol: NPHS2",
        "   Full name: NPHS2 stomatin family member, podocin",
        "   Aliases: SRN1, PDCN",
        "   HG38 coordinates: 1:179550539-179575952",
        "   HG19 coordinates: 1:179519674-179545087",
        # "6. HERC2 (HGNC:4868)"
        "   Symbol: HERC2",
        "   Full name: HECT and RLD domain containing E3 ubiquitin protein ligase 2",
        "   Aliases: jdf2, p528, D15F37S1",
        "   HG38 coordinates: 15:28111040-28322179",
        "   HG19 coordinates: 15:28356186-28567298",
    ]

    lines = result.output.splitlines()

    for e in expected:
        assert e in lines


@pytest.mark.parametrize(
    "invalid_pmid",
    [
        "invalid",
        "0",
        "01234",
        "123456789",  # Too long
        "",
    ],
)
def test_cli_rejects_invalid_pmid(cli_runner, invalid_pmid):
    """Test CLI rejects invalid PMID formats."""
    result = cli_runner.invoke(main, ["--pmid", invalid_pmid])
    assert result.exit_code != 0
    assert "Invalid PMID" in result.output
