"""
Test CLI functionality
"""

import pytest
from click.testing import CliRunner

from pmid_parser.main import main


@pytest.fixture
def cli_runner():
    """Provide a Click CLI test runner."""
    return CliRunner()


def test_cli_requires_pmid(cli_runner):
    """Test that CLI requires PMID argument."""
    result = cli_runner.invoke(main, [])
    assert result.exit_code != 0
    assert "Missing option '--pmid'" in result.output


def test_cli_works(cli_runner):
    """Test CLI accepts PMID and doesn't crash."""
    result = cli_runner.invoke(main, ["--pmid", "38790019"])
    assert result.exit_code == 0
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
    print(result.output.splitlines()[0])
    assert result.output.splitlines()[0].strip() == "Processing PMID: 38790019"


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
