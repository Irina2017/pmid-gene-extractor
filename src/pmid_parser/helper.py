"""
Helper functions
"""

import csv
from typing import List

from .gene_extractor import ValidatedGeneRecord


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


def write_genes_to_csv(
    validated_genes: List[ValidatedGeneRecord], filename: str
) -> None:
    """Write validated genes to CSV file."""

    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)

        # Write header
        writer.writerow(
            [
                "HGNC ID",
                "HGNC Gene Name",
                "Gene Aliases",
                "Hg38 genomic coordinates",
                "Hg19 genomic coordinates",
                "Disease",
            ]
        )

        # Write gene data
        for _i, gene in enumerate(validated_genes, 1):
            writer.writerow(
                [
                    gene.hgnc_id,
                    gene.hgnc_symbol,
                    gene.gene_aliases,
                    gene.hg38_coordinates,
                    gene.hg19_coordinates,
                    " | ".join(gene.diseases) if gene.diseases else "",
                ]
            )
