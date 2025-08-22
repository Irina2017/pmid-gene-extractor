"""
Helper functions
"""


def is_valid_pmid(pmid: str) -> bool:
    """
    Function to check whether pmid string is valid PMID identifier
    Args:
        pmid: string to validate

    Returns
        True if valid format; otherwise, False

    NOTE: valid format is
           -  1 to 8 digits
           -  no leading 0's
           -  if treated as number, then positive (not a 0)
    """

    # not empty string
    if not pmid:
        return False

    # check if all characters are digits
    if not pmid.isdigit():
        return False

    # no leading 0 and cannot be 0 itself
    if pmid.startswith("0"):
        return False

    # length is between 1 and 8 digits
    return len(pmid) <= 8
