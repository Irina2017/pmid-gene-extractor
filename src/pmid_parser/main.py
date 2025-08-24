"""
Command-line interface for PMID parser
"""

from typing import Optional

import click

from .gene_extractor import extract_and_validate_genes_from_paper
from .helper import is_valid_pmid, write_genes_to_csv
from .pmid_reader import read_paper_content


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

    # read paper by PMID
    try:
        paper_content = read_paper_content(pmid)

        # Display results
        click.echo(f"Title: {paper_content['metadata']['title']}")
        click.echo(f"Abstract: {len(paper_content['metadata']['abstract'])} characters")
        click.echo(f"PMC ID: {paper_content['metadata']['pmc_id']}")
        # click.echo(f"MeSH terms: {len(paper_content['metadata']['mesh_terms'])} found")
        # click.echo(f"MeSH terms: {paper_content['metadata']['mesh_terms']} found")
        # click.echo(f"Sections found: {list(paper_content['sections'].keys())}")

        # section previews
        for section_name, content in paper_content["sections"].items():
            if content:
                click.echo(f"{section_name}: {len(content)} paragraphs")

        # After paper content extraction:
        gene_results = extract_and_validate_genes_from_paper(paper_content)

        click.echo("VALIDATED GENE RESULTS:")
        click.echo(f"Found {len(gene_results['validated_genes'])} validated genes")

        for i, gene in enumerate(gene_results["validated_genes"], 1):
            click.echo(f"{i}. {gene.hgnc_symbol} ({gene.hgnc_id})")
            click.echo(f"   Symbol: {gene.hgnc_symbol}")
            click.echo(f"   Full name: {gene.full_name}")
            click.echo(f"   Aliases: {gene.gene_aliases}")
            click.echo(f"   HG38 coordinates: {gene.hg38_coordinates}")
            click.echo(f"   HG19 coordinates: {gene.hg19_coordinates}")
            click.echo(f"   Mentioned in: {', '.join(list(set(gene.sources)))}")
            click.echo(f"   Diseases: {gene.diseases}")
            click.echo()

        if output:
            write_genes_to_csv(gene_results["validated_genes"], output)

        click.echo("Paper content extracted successfully!")

    except Exception as e:
        raise click.ClickException(str(e)) from e
