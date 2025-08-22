"""
Command-line interface for PMID parser
"""

from typing import Optional

import click

from .helper import is_valid_pmid


@click.command()
@click.option(
    "--pmid", required=True, help="PubMed ID to parse (e.g., 38790019 or 'PMID: 38') "
)
@click.option("--output", "-o", help="Path to output file (default:stdout)")
@click.version_option()
def main(pmid: str, output: Optional[str]) -> None:
    """
    Extract genes, co-mentioned diseases from scientific text; added extra metadata

    Args:

         - pmid:   PubMed ID to parse, e.g. 38790019, or PMID: 38790019.

         - output: path to file to store results (optional)

    Returns:

         None

    """

    # clean and validate pmid
    pmid = pmid.strip()
    if pmid.startswith("PMID:"):
        pmid = pmid[5:].strip()
    if not is_valid_pmid(pmid):
        raise click.BadParameter(f"Invalid PMID: {pmid}")

    click.echo(f"Processing PMID: {pmid}")

    # read PMID

    if output:
        click.echo(f"Output will be saved to: {output}")
    else:
        click.echo("Output will be printed to stdout")
    click.echo("CLI scaffold working!")
