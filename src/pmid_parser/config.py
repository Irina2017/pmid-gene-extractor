"""
Configuration constants for PMID parser.
"""

# API endpoints
PUBMED_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PMC_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


# API settings
API_TIMEOUT = 30
# Without API key, NCBI's E-utilites allow 3 requests per seconds =>
#   need only 2 API calls - no need
#   In batch/production, would add an API key + backoff and consider lowering the delay.
RATE_LIMIT_DELAY_SECONDS = 0.5  # needed for pytests

# ENSEMBL settings
ENSEMBL_REST = "https://rest.ensembl.org"
ENSEMBL_GRCH37_REST = "https://grch37.rest.ensembl.org"
ENSEMBL_HEADERS = {"Accept": "application/json", "User-Agent": "pmid-gene-parser/0.1"}


# REFERENCES
HGNC = "references/HGNC_download_custom_Aug2025.txt"
