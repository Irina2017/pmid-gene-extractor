"""
Configuration constants for PMID parser.
"""

# API endpoints
PUBMED_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PMC_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


# API settings
API_TIMEOUT = 30
# Without API key, NCBI's E-utilites allow 3 requests per seconds =>
#   use 0.5s delay (~2 requests per second) as a safe buffer for one PMID assessment.
#   In batch/production, would add an API key + backoff and consider lowering the delay.
RATE_LIMIT_DELAY_SECONDS = 0.5
