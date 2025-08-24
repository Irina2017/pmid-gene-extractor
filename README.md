# pmid-gene-extractor

# PMID Gene Extractor

A Python tool for extracting gene mentions and disease associations from scientific publications using PubMed IDs.

## Goal

Extract genes in HGNC format and their disease associations from biomedical literature by:
- Fetching full-text articles from PubMed/PMC when available
- Parsing structured sections (title, abstract, methods, results, discussion)
- Identifying gene mentions and validating against HGNC database
- Mapping gene-disease relationships with provenance tracking

## Installation

This project uses Poetry for dependency management:

```bash
# Clone the repository
git clone https://github.com/Irina2017/pmid-gene-extractor.git
cd pmid-gene-extractor

# Install dependencies
poetry install

# Activate the environment
poetry shell
```

## Usage


```bash
# Extract genes from a paper
poetry run pmid-gene-parser --pmid 38790019

# Save results to file
poetry run pmid-gene-parser --pmid 38790019 --output results.csv

#help
poetry run pmid-gene-parser --help
```

## Project Structure

src/pmid_parser/main.py - Main entry point and CLI interface
src/pmid_parser/pmid_reader.py - PubMed/PMC API integration and text extraction
src/pmid_parser/gene_extractor.py - Extraction of gene candidates, validation against HGNC DB, disease extraction
src/pmid_parser/helpers.py - Utility functions for validation
src/pmid_parser/config.py - Configuration constants


## Documentation
Detailed documentation is available in the docstrings of each module. The main entry point is main.py.

## Requirements
Python 3.12+
Internet access for PubMed/PMC APIs
Dependencies managed via pyproject.toml
