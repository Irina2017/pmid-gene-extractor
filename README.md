# PMID Gene Extractor

A Python tool for extracting gene mentions and disease associations from scientific publications using PubMed IDs.

## Goal

Extract genes in HGNC format and their disease associations from biomedical literature by:
- Fetching full-text articles from PubMed/PMC when available
- Parsing structured sections (title, abstract, methods, results, discussion)
- Identifying gene mentions and validating against HGNC database
- Mapping gene-disease relationships

## Documentation
Detailed documentation is available in the docstrings of each module. The main entry point is main.py.

## Requirements
 -- Python 3.12+

 -- Internet access for PubMed/PMC APIs

 -- Dependencies managed via pyproject.toml

## Installation

This project uses Poetry for dependency management:

```bash
# Clone the repository
git clone https://github.com/Irina2017/pmid-gene-extractor.git
cd pmid-gene-extractor

# Install poetry - either through curl or pip
#(in later case - comment curl command and uncomment pip command)
curl -sSL https://install.python-poetry.org | python3 -
#pip install poetry


# Install dependencies
poetry install

```

### NOTE

 If during installation poetry complains about Python version, then:
```bash

# Install Homebrew if needed: https://brew.sh/

# Install pyenv if not installed
brew install pyenv

# Install specific Python version
pyenv install 3.12.7

# Set as global default
pyenv global 3.12.7

# Install dependencies
poetry install

```


## Usage

```bash
# Extract genes from a paper
# Save results to file - produces csv file with expected fields - HGNC id, HGNC gene symbol,
#aliases, HG38 and HG19 coordinates and extracted disease
poetry run pmid-gene-parser --pmid 38790019 --output results.csv


# Without output option, progam will output some statistics about the genes
# but in different format than above
poetry run pmid-gene-parser --pmid 38790019

#help
poetry run pmid-gene-parser --help
```

### NOTE

If you see:

```bash
urllib3 v2 only supports OpenSSL 1.1.1+, currently the 'ssl' module is compiled with 'LibreSSL ...'
```
your Python is linked against LibreSSL. Use a Python build with OpenSSL if you can.
That warning does not prevent program to run successfully though.

## Project Structure

-- src/pmid_parser/main.py - Main entry point and CLI interface

-- src/pmid_parser/pmid_reader.py - PubMed/PMC API integration and text extraction

-- src/pmid_parser/gene_extractor.py - Extraction of gene candidates, validation against HGNC DB, disease extraction

-- src/pmid_parser/helpers.py - Utility functions for validation

-- src/pmid_parser/config.py - Configuration constants


