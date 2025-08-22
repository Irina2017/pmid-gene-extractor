"""
Test helper functions.
"""

import pytest

from pmid_parser.helper import is_valid_pmid


@pytest.mark.parametrize(
    "pmid, expected",
    [
        # valid
        ("1", True),
        ("38790019", True),
        ("12345678", True),
        # Empty string
        ("", False),
        # Non-digit characters
        ("12a34", False),
        ("-1234", False),
        # Leading zeros
        ("0", False),
        ("00000001", False),
        # Too long (>8 digits)
        ("123456789", False),
    ],
)
def test_is_valid_pmid(pmid: str, expected: bool) -> None:
    """
    Tests for is_valid_pmid() from pmid_extractor.helpers
    """
    assert is_valid_pmid(pmid) == expected
