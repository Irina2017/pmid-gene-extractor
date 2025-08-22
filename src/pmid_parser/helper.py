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
        print(f"ERROR: not all digits in PMID - {pmid} ")
        return False

    # no leading 0 and cannot be 0 itself
    if pmid.startswith("0"):
        print(f"ERROR: PMID has leading 0 - {pmid} ")
        return False

    # length is between 1 and 8 digits
    if len(pmid) > 8:
        print(f"ERROR: PMID has more than 8 digits - {pmid} ")
        return False

    return True
