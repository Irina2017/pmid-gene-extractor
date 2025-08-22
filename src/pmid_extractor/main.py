"""
Command-line interface for PMID parser
"""

from typing import Optional

import click


@click.command()
@click.option("--pmid", required=True, help="PubMed ID to process (e.g., 38790019) ")
@click.option("--output", "-o", help="Path to output file (default:stdout)")
@click.version_option()
def main(pmid: str, output: Optional[str]) -> None:
    "Extract genes, co-mentioned diseases from scientific text; added extra metadata"
    click.echo(f"Processing PMID: {pmid}")
    if output:
        click.echo(f"Output will be saved to: {output}")
    else:
        click.echo("Output will be printed to stdout")
    click.echo("CLI scaffold working!")
