"""Simple CLI test: (simplified due time-constrants)"""

from click.testing import CliRunner

from pmid_extractor.main import main


def test_cli_works():
    """Test CLI accepts PMID and doesn't crash."""
    runner = CliRunner()
    result = runner.invoke(main, ["--pmid", "38790019"])
    assert result.exit_code == 0
