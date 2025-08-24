"""
gene_extractor: module to extract and validate gene mentions from available sections/metadata (assessment MVP).

Core Goal: Extract ≥5 HGNC-validated gene symbols from paper content for genetic research analysis.

- Extract gene candidates using regex backup
- Validate candidates against HGNC database
- Return validated genes with source section reference (title, abstract, results, etc.)

Gene Extraction Pipeline:

NOTE:  spaCy models (en_ner_bionlp13cg_md, en_ner_jnlpba_md) had poor performance on the test PMID: 38790019. => removed

Current decision: For a short paper (<10 pages, single PMID), process all sections, combine findings, validate against HGNC,
and warn if fewer than 5 unique genes are found. If more than 5 are validated, optionally filter or score by section importance
or mention frequency if needed.

Alternative(longer texts): process sections in priority order, validating after each. Once a gene is confirmed, scan all sections
for its HGNC ID, symbol, aliases, and full name to capture mentions. Continue regex on new sections only until at least 5 genes are
validated; after that, restrict to exact matches of confirmed genes.

1. Regex pattern matching  - several regex to catch either HGNC id or gene symbols in HGNC format.
    Currently, it does not catch cases like BRCA1/2 or "beta-carotene oxygenase 2":
       for MVP, focus is on gene id and gene symbols as they seen in HGNC database, since
         at least 5 genes (not all) were requested.
    Regex looks for genes only ( no pseudogenes, no long non-coding RNAs or any RNAs))

    Several regex may return the same gene, so duplicates are removed later - fine for one PMID.
    Point of tuning later.

Sections are processed in the next order:
- Title:  main genes studied
- Abstract: key findings summary
- Intro: background details
- Results: experimental findings (cases are included here)
- Discussion: interpretation, related genes
- Methods: technical notes
- Conclusions: summary statements

NOTE: In future, abbreviations from the paper’s abbreviation section should be added to STOP_KEYWORDS to reduce false positives.
For example, MGP may appear both as an abbreviation and as a valid gene symbol. Extra context analysis would be needed
to distinguish such cases. In this MVP, MGP is not scored because ≥5 valid genes are already found without it, but it would still
be processed if referenced by its HGNC ID.

2.  All candidates are verified against HGNC REST API/DB

Uses pre-downloaded HGNC database (references/HGNC_download_custom_Aug2025.txt) for gene validation instead of live API calls.

Database approach benefits vs API calls:
- Performance: Instant validation (milliseconds) vs 10+ seconds for API calls
- Reliability: Zero network dependency, deterministic results
- Simplicity: Dictionary lookups vs rate limiting, error handling, timeouts


Database limitations (acceptable for MVP):
- Static data: No real-time updates (HGNC changes are infrequent)
- Storage: ~5-10MB in git repo vs zero for API approach
- Setup: Requires data file download

Database source: HUGO Gene Nomenclature Committee (HGNC) - official human gene naming authority.
File location: references/HGNC_download_custom_Aug2025.txt (included in repository for assessment convenience)

a) validate all with HGNC id
b) validate by potential gene name vs approved symbol
if more than 5 genes are validated, skip validating the rest vs aliases,previous names etc (outside of scope MVP)


3. NOTE: Once all information is received from HGNC for confirmed entities, another search in text can be done
with gene names like "beta-carotene oxygenase 2" for those confirmed to ensure all mentions are captured - did not benefit in this case =>
REMOVED.

4. -> 3. Disease/Disorder Association Extraction:
   After gene validation, sentences containing validated genes are analyzed for disease associations
   using token-based pattern matching.
   (1) Gene-related patterns ( "RRAGD-related hypomagnesemia") and
   (2) Association keywords within 20-token window .
   Uses parenthesis-aware parsing to extract disease phrases from gene mentions to  sentence delimiters.
   Filters sentences to keep only those with clear disease associations, providing gene-disease relationship context for each validated gene.

Error handling (MVP):
- No genes found -> return empty list with warning
- <5 genes found -> log warning, return all found genes

TODO: If time permits - move some of the logic to separate modules:
        regex_patterns.py      # All regex extraction functions
        hgnc_validation.py     # HGNC database functions
        disease_extraction.py  # Disease-related functions
        ensembl_lookup.py      # Ensembl coordinate functions


1. MAIN GENE EXTRACTION
extract_and_validate_genes_from_paper() - Main entry point
extract_gene_candidates()
_process_section_with_tracking()
_extract_with_regex()
_consolidate_sentence_findings()
_remove_stop_words_from_sentence()

2. HGNC VALIDATION
load_hgnc_database()
validate_gene_candidates
validate_hgnc_ids()
validate_gene_symbols()
_merge_records()
_create_validated_gene_records()
partition_matches_by_hgnc()
partition_matches_by_symbol()
merge_temp_to_lists()

3. DISEASE EXTRACTION
extract_disease_associations
analyze_disease_sentences_token_based()
extract_diseases_from_sentence()
analyze_disease_sentences_list()
_extract_until_delimiter()
_check_gene_related_pattern()
_check_association_pattern()

4. ENSEMBL LOOKUP
_lookup_ensembl_coordinates()
_lookup_gene()
_http_get()
"""

import re
from dataclasses import dataclass

# import time
from typing import Any, Dict, List, Optional, Set, Tuple

import requests

from .config import (
    API_TIMEOUT,
    ENSEMBL_GRCH37_REST,
    ENSEMBL_HEADERS,
    ENSEMBL_REST,
    HGNC,
)

# stop keywords to ignore (common false positives)
# TODO - later pull abbreviations from the text and add them here to exclude as well
STOP_KEYWORDS: Set[str] = {
    # Methods / tech/ db /regulations
    "DNA",
    "RNA",
    "PCR",
    "NGS",
    "WGS",
    "WES",
    "CNV",
    "SNP",
    "SNV",
    "IHC",
    "EDTA",
    "PBS",
    "ATP",
    "ADP",
    "NAD",
    "NADH",
    "CI",
    "SD",
    "BAM",
    "CRAM",
    "VCF",
    "FASTQ",
    "FASTA",
    "GTF",
    "BED",
    "TSV",
    "CSV",
    "JSON",
    "VUS",
    "LP",
    "GS",
    "ES",
    "EGBP",
    "MIM",
    "HIPAA"
    # Figures / refs
    "FIGURE",
    "TABLE",
    "FIG",
    "TAB",
    "ETAL",
    # Controls / genotypes
    "WT",
    "KO",
    # Organizations / regions
    "USA",
    "US",
    "UK",
    "EU",
    "WHO",
    "FDA",
    "NIH",
    "CDC",
    "IRB",
    "CLIA",
    "CAP",
    # Misc tokens often misfired as genes
    "MT",
    "HLA",
    # Stopwords (shouty case from PDFs)
    "AND",
    "OR",
    "NOT",
    "THE",
    "ALL",
    "ANY",
    "BUT",
    "FOR",
    "TO",
    "IN",
    "ON",
    # current abbreviation in the paper - later pull from the paper for automatic adjustment
    # if another paper is processed
    "MGP",
    "FSGS",
}


def load_hgnc_database():
    """Load HGNC CSV into dictionaries - fastest lookup."""
    hgnc_by_id = {}
    hgnc_by_symbol = {}

    with open(HGNC, "r") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            fields = line.strip().split("\t")
            row = dict(zip(header, fields, strict=False))

            hgnc_by_id[row["hgnc_id"]] = row
            hgnc_by_symbol[row["approved_symbol"]] = row

    return hgnc_by_id, hgnc_by_symbol


# Load once at module startup
HGNC_BY_ID, HGNC_BY_SYMBOL = load_hgnc_database()


@dataclass
class GeneMatch:
    """Structure for regex-based gene findings."""

    gene_name: str
    gene_id: str = ""  # HGNC ID if found
    section: str = ""
    element_index: int = 0
    sentence_index: int = 0
    sentence: str = ""
    full_element: str = ""  # Complete element text
    method: str = "regex"


@dataclass
class ValidatedGeneRecord:
    """Complete gene record with HGNC validation and mentions."""

    hgnc_id: str
    hgnc_symbol: str
    gene_aliases: List[str]
    full_name: str
    hg38_coordinates: str = ""
    hg19_coordinates: str = ""
    ensembl_id: str = ""
    omim_id: str = ""
    diseases: List[str] = None  # Empty for now
    mentions: List[str] = None  # All sentences where mentioned (deduplicated)
    sources: List[str] = None  # Which sections mentioned
    full_paragraph: List[str] = None  # which paragraph

    def __post_init__(self):
        if self.diseases is None:
            self.diseases = []
        if self.mentions is None:
            self.mentions = []
        if self.sources is None:
            self.sources = []


def extract_and_validate_genes_from_paper(paper_content: Dict) -> Dict:
    """
    Extract and validate genes from paper with HGNC database integration.

    Args:
        paper_content: Dict with 'metadata' and 'sections' keys

    Returns:
        Dict with validated gene records and raw findings
    """
    results = {
        "regex_findings": [],  # records that fail validation for manual review
        "validated_genes": [],  # ValidatedGeneRecord objects
        "timing": {},
        "section_summaries": {},
    }

    # Step 1: Extract gene candidates from all sections
    found_hgnc_ids = extract_gene_candidates(paper_content, results)

    # Step 2: Validate candidates against HGNC database
    validated_records = validate_gene_candidates(found_hgnc_ids, results)

    # Step 3: Extract disease associations from validated genes
    extract_disease_associations(validated_records)

    results["validated_genes"] = validated_records
    return results


def extract_gene_candidates(paper_content: Dict, results: Dict) -> Set[str]:
    """Extract gene candidates from all paper sections."""

    found_hgnc_ids = set()
    section_order = [
        "title",
        "abstract",
        "intro",
        "results",
        "discussion",
        "methods",
        "conclusion",
    ]

    print("=" * 80)
    print("CANDIDATES EXTRACTION FROM PAPER")

    # Process metadata sections first (title, abstract)
    for section_name in ["title", "abstract"]:
        if section_name in paper_content["metadata"]:
            content = [paper_content["metadata"][section_name]]
            section_findings = _process_section_with_tracking(
                section_name, content, found_hgnc_ids
            )
            results["regex_findings"].extend(section_findings)

    # Process body sections
    for section_name in section_order[2:]:
        if section_name in paper_content.get("sections", {}):
            content = paper_content["sections"][section_name]
            section_findings = _process_section_with_tracking(
                section_name, content, found_hgnc_ids
            )
            results["regex_findings"].extend(section_findings)

    return found_hgnc_ids


def validate_gene_candidates(
    found_hgnc_ids: Set[str], results: Dict
) -> List[ValidatedGeneRecord]:
    """Validate gene candidates against HGNC database."""

    print("=" * 80)
    print("VALIDATION AGAINST HGNC")

    # a) Validate HGNC IDs
    validated_hgnc_records = validate_hgnc_ids(found_hgnc_ids)
    (validated_hgnc_records, results["regex_findings"]) = _merge_records(
        validated_hgnc_records, results["regex_findings"]
    )

    # b) Validate remaining gene symbols
    remaining_symbols = {f.gene_name for f in results["regex_findings"] if f.gene_name}
    validated_symbol_records = validate_gene_symbols(remaining_symbols)
    (validated_symbol_records, results["regex_findings"]) = _merge_records(
        validated_symbol_records, results["regex_findings"]
    )

    # Check if we have enough genes
    total_validated = len(validated_hgnc_records) + len(validated_symbol_records)
    if total_validated < 5:
        print("WARNING: Less than 5 genes were validated from the paper by this MVP.")
        remaining_symbols = {
            f.gene_name for f in results["regex_findings"] if f.gene_name
        }
        if remaining_symbols:
            print(
                "Future development - to verify remaining potential gene symbols vs aliases and previous names."
            )
            print(f"Left potential gene symbols to check: {remaining_symbols}")
        else:
            print("No potential gene symbols left to check.")

    # Create validated gene records
    validated_records = _create_validated_gene_records(
        validated_hgnc_records | validated_symbol_records
    )

    return validated_records


def extract_disease_associations(validated_records: List[ValidatedGeneRecord]) -> None:
    """Extract disease associations from validated gene mentions."""

    print("=" * 80)
    print("TOKEN-BASED DISEASE EXTRACTION")

    for gene_record in validated_records:
        analyze_disease_sentences_token_based(gene_record)


def _http_get(url: str) -> Any:
    r = requests.get(url, headers=ENSEMBL_HEADERS, timeout=API_TIMEOUT)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


def _lookup_gene(server: str, ensembl_id: str) -> Optional[Dict[str, Any]]:
    # GET /lookup/id/:id  (returns seq_region_name/start/end/strand/external_name/assembly_name/object_type)
    url = f"{server}/lookup/id/{ensembl_id}?expand=0"
    data = _http_get(url)
    if not data or data.get("object_type") != "Gene":
        return None
    return data


def _merge_records(
    records: Dict[str, Dict], findings: List[GeneMatch]
) -> Tuple[Dict[str, Dict], List[GeneMatch]]:
    """
    Merge mentions etc per HGNC id/or HGNS symbol
    """
    for v in records.keys():
        temp, left_to_process = partition_matches_by_hgnc(
            findings, records[v]["hgnc_id"], records[v]["approved_symbol"]
        )
        (sentences, sources, full) = merge_temp_to_lists(temp)
        pruned_sentences = analyze_disease_sentences_list(
            sentences, records[v]["hgnc_id"], records[v]["approved_symbol"]
        )
        records[v][
            "sentences"
        ] = pruned_sentences  # need to keep only those that have "diseases"
        records[v]["sources"] = sources
        records[v][
            "full_elements"
        ] = full  # currently keep as level of truth, may remove later

        findings = (
            left_to_process  # keep only unvalidated in  results['regex_findings']
        )
    return (records, findings)


def _extract_until_delimiter(
    tokens: List[str], start_idx: int, trigger_idx: int
) -> List[str]:
    """
    Extract tokens from start_idx with parenthesis-aware stopping logic (like Polish notation idea)

    Args:
        start_idx: index where gene mention starts
        trigger_idx: index where 'related' or association keyword was found

    Logic:
    1. From start_idx to trigger_idx: count parentheses (+ for '(', - for ')')
    2. After trigger_idx: stop at ')' only if parenthesis count goes negative
    3. Always stop at '.'
    """

    extracted_tokens = []
    paren_count = 0

    for i in range(start_idx, len(tokens)):
        token = tokens[i]
        extracted_tokens.append(token)

        # count parentheses in current token
        paren_count += token.count("(")
        paren_count -= token.count(")")

        # after passing trigger index (related/associated keyword)
        if i > trigger_idx:
            # Stop at closing parenthesis only if count goes negative
            if token.endswith(")") and paren_count == -1:
                return extracted_tokens

        # Always stop at dot
        if token.endswith(
            "."
        ):  # problem if dot is not end of sentence, but part of NM_014625.3 - for later
            return extracted_tokens

    return extracted_tokens


def _check_gene_related_pattern(
    tokens: List[str], gene_idx: int, gene_symbol: str
) -> str:
    """
    Check for gene-related pattern within next 3 tokens.
    Cases like RRAGD-related
    """

    for i in range(gene_idx, min(gene_idx + 4, len(tokens))):
        if "related" in tokens[i].lower():
            disease_tokens = _extract_until_delimiter(tokens, gene_idx, i)
            disease_string = " ".join(disease_tokens).strip()
            return disease_string

    return ""


def _check_association_pattern(
    tokens: List[str], gene_idx: int, disease_keywords: List[str]
) -> str:
    """Check for association keywords within window of 20 tokens."""

    window_start = max(0, gene_idx - 10)
    window_end = min(len(tokens), gene_idx + 11)

    # Look for disease keywords in window
    for i in range(window_start, window_end):
        token_clean = tokens[i].lower().strip(".,()[]{}:;!?")
        if token_clean in disease_keywords:
            disease_tokens = _extract_until_delimiter(tokens, gene_idx, i)
            disease_string = " ".join(disease_tokens).strip()
            return disease_string
    return ""


def extract_diseases_from_sentence(
    sentence: str, gene_symbol: str, gene_id: str
) -> List[str]:
    """
    Extract diseases from sentence using token-based analysis.

    Args:
        sentence: sentence to analyze
        gene_symbol: gene symbol
        gene_id:    gene HGNC id

    Returns:
        List of potential disease strings
    """

    disease_keywords = [
        "disease",
        "syndrome",
        "disorder",
        "symptoms",
        "symptom",
        "association",
        "associated",
        "related",
        "cause",
        "caused",
        "condition",
    ]

    tokens = sentence.split()
    diseases_found = []

    # Find all gene mention indices
    gene_indices = []
    for i, token in enumerate(tokens):
        clean_token = token.lower()
        if (
            clean_token == gene_symbol.lower()
            or clean_token == gene_id.lower()
            or gene_symbol.lower() in clean_token
            or gene_id.lower() in clean_token  # like "RRAGD-related"
        ):
            gene_indices.append(i)

    # For each gene mention, look for diseases
    for gene_idx in gene_indices:
        # CASE 1: Gene + "related" within next 3 tokens
        case1_disease = _check_gene_related_pattern(tokens, gene_idx, gene_symbol)
        if case1_disease:
            diseases_found.append(case1_disease)
            # print(f"   CASE 1 - Gene-related: '{case1_disease}'")
            continue

        # CASE 2: Association keywords within window of 20 tokens
        case2_disease = _check_association_pattern(tokens, gene_idx, disease_keywords)
        if case2_disease:
            diseases_found.append(case2_disease)
            # print(f"   CASE 2 - Association: '{case2_disease}'")

    return diseases_found


def analyze_disease_sentences_token_based(validated_gene: ValidatedGeneRecord) -> None:
    """
    Analyze sentences using token-based approach.
    See if associated/association or <gene-name>-related or disease keywords are within small window
    If yes - pick up those substrings

    """

    gene_symbol = validated_gene.hgnc_symbol
    gene_id = validated_gene.hgnc_id

    all_diseases = []
    kept_sentences = []

    # print("="*60)
    # print(f"{gene_symbol}")

    for _i, sentence in enumerate(validated_gene.mentions, 1):
        diseases = extract_diseases_from_sentence(sentence, gene_symbol, gene_id)

        if diseases:
            all_diseases.extend(diseases)
            kept_sentences.append(sentence)

    # Update gene record
    validated_gene.mentions = kept_sentences
    validated_gene.diseases = list(set(all_diseases))  # Remove duplicates

    # print(f"SUMMARY for {gene_symbol}:")
    # print(f"   Diseases extracted: {validated_gene.diseases}")
    # print(f"   Sentences: {len(validated_gene.mentions)}/{len(validated_gene.mentions)} kept")


def analyze_disease_sentences_list(
    sentences: List[str], hgnc_id: str, gene_symbol: str
) -> List[str]:
    """Analyze sentences for disease associations and keep only relevant ones"""

    original_count = len(sentences)

    print(f"ANALYZING {gene_symbol} - Original mentions: {original_count}")

    case1_sentences = []  # gene-related disease patterns
    case2_sentences = []  # disease keyword patterns
    removal_sentences = []  # no disease association

    # Case 1: Look for <gene_name>-related patterns (FIXED - capture only the disease part)
    gene_related_pattern = (
        rf"\b{re.escape(gene_symbol.lower())}-related\s+([^,)]*?)(?:[,)]|$)"
    )

    # Case 2: Disease/medical keywords (FIXED - more specific)
    disease_keywords = [
        "disease",
        "syndrome",
        "disorder",
        "symptoms",
        "symptom",
        "association",
        "associated with",
        "related to",
        "cause",
        "caused by",
        "condition",
        "associated",
    ]

    for _i, sentence in enumerate(sentences, 1):
        sentence_clean = sentence.strip()
        sentence_lower = sentence_clean.lower()

        # is gene mentioned in this sentence
        gene_mentioned = (
            gene_symbol.lower() in sentence_lower or hgnc_id.lower() in sentence_lower
        )

        if not gene_mentioned:
            continue

        # check for gene-related disease patterns
        case1_match = re.search(gene_related_pattern, sentence_lower)
        if case1_match:
            # extracted_disease = case1_match.group(1).strip()
            case1_sentences.append(sentence_clean)
            continue

        # check for disease keywords in same sentence as gene
        has_disease_keyword = any(
            keyword in sentence_lower for keyword in disease_keywords
        )
        if has_disease_keyword:
            # found_keywords = [kw for kw in disease_keywords if kw in sentence_lower]
            case2_sentences.append(sentence_clean)
            continue

        # mark for removal - neither of cases
        removal_sentences.append(sentence_clean)

    # update gene record with filtered sentences
    filtered_sentences = case1_sentences + case2_sentences
    final_count = len(filtered_sentences)

    print(f"   Removed (no disease): {len(removal_sentences)} sentences")
    print(f"   Total: {original_count} -> {final_count} sentences")
    print("=" * 60)

    return filtered_sentences


def _lookup_ensembl_coordinates(ensembl_id: str, url: str) -> str:
    """
    Function to lookup Ensembl coordinates by EnsemblID

    Args:
         ensembl_id: string to query Ensemble ID API
         url: must be either ENSEMBL_GRCH37_REST or ENSEMBL_REST values

    Returns:
         chr:start-end if Ensembl returns record
         Otherwise - give an error message and returns ''
    """

    chrom = start_pos = end_pos = ""
    try:
        genome_coordinates = _lookup_gene(url, ensembl_id)
        if not genome_coordinates:
            raise ValueError(
                f"Ensembl lookup  failed for ensembl_id : {ensembl_id}, at {url}"
            )

        chrom = genome_coordinates["seq_region_name"]
        start_pos = genome_coordinates["start"]
        end_pos = genome_coordinates["end"]
        # strand = genome_coordinates["strand"]
    except ValueError as e:
        print(f"Warning: {e}")

    return f"{chrom}:{start_pos}-{end_pos}" if genome_coordinates else ""


def _create_validated_gene_records(
    records_dict: Dict[str, Dict],
) -> List[ValidatedGeneRecord]:
    """Create ValidatedGeneRecord objects from records dictionary."""
    validated_genes = []

    for rec in records_dict.values():
        gene = ValidatedGeneRecord(
            hgnc_id=rec["hgnc_id"],
            hgnc_symbol=rec["approved_symbol"],
            gene_aliases=rec.get("alias_symbols", []),
            full_name=rec.get("approved_name", ""),
            ensembl_id=rec.get("ensembl_gene_id", ""),
            hg38_coordinates=rec.get("genomic_pos_hg38", ""),
            hg19_coordinates=rec.get("genomic_pos_hg19", ""),
            diseases=[],  # Placeholder, no disease extraction in MVP
            omim_id=rec.get("omim_id", ""),
            mentions=rec.get("sentences", []),
            sources=rec.get("sources", []),
            full_paragraph=rec.get("full_elements", []),
        )

        gene.hg38_coordinates = _lookup_ensembl_coordinates(
            gene.ensembl_id, ENSEMBL_REST
        )
        gene.hg19_coordinates = _lookup_ensembl_coordinates(
            gene.ensembl_id, ENSEMBL_GRCH37_REST
        )

        validated_genes.append(gene)

    return validated_genes


def merge_temp_to_lists(
    temp: List[GeneMatch],
) -> Tuple[List[str], List[str], List[str]]:
    """
    From a temp list of GeneMatch:
      - keep only unique (sentence, source) pairs
      - return (sentences, sources, full) as lists
    Source uniqueness uses (section, element_index), but 'sources' returned
    is just the section name per kept item (simple list, no compression).
    """
    # 1) de-dup by (sentence, source)
    seen_pairs = set()
    deduped = []  # store (sentence_norm, section_name, full_text)
    for m in temp:
        sent_norm = " ".join((m.sentence or "").split()).strip()
        if not sent_norm:
            continue
        section = (m.section or "").strip()
        src_key = (section, m.element_index)  # source identity
        key = (sent_norm, src_key)  # (sentence, source)
        if key in seen_pairs:
            continue
        seen_pairs.add(key)
        deduped.append((sent_norm, section, m.full_element or ""))

    # 2) build simple lists (no extra uniq/sorting)
    sentences = [s for s, _, _ in deduped]
    sources = [sec for _, sec, _ in deduped]  # simple list of section names
    full = [f for _, _, f in deduped if f]

    return (sentences, sources, full)


def partition_matches_by_hgnc(
    results: List[GeneMatch],
    hgnc_id: str,
    symbol: str,
) -> Tuple[List[GeneMatch], List[GeneMatch]]:
    """
    Split results['regex_findings'] into:
      - temp: records that match any of:
          (a) same gene_id == hgnc_id AND same gene_name == symbol
          (b) same gene_id == hgnc_id
          (c) same gene_name == symbol AND empty gene_id
      - left_to_process: all other records

    Returns (temp, left_to_process).
    """
    # Normalize targets
    """
    temp: records matching any of:
      (a) gene_id == hgnc_id AND gene_name == symbol
      (b) gene_id == hgnc_id
      (c) gene_name == symbol AND gene_id == ""
    left_to_process: all others
    """
    # normalize inputs
    target_id = (hgnc_id or "").strip().upper()
    if target_id and not target_id.startswith("HGNC:"):
        target_id = f"HGNC:{target_id}"
    target_symbol = (symbol or "").strip().upper()

    temp: List[GeneMatch] = []
    left: List[GeneMatch] = []

    for gm in results:
        gm_id = (gm.gene_id or "").strip().upper()
        gm_sym = (gm.gene_name or "").strip().upper()

        cond_a = gm_id == target_id and gm_sym == target_symbol
        cond_b = gm_id == target_id
        cond_c = gm_sym == target_symbol and gm_id == ""

        (temp if (cond_a or cond_b or cond_c) else left).append(gm)

    return temp, left


def partition_matches_by_symbol(
    results: List[GeneMatch],
    symbol: str,
) -> Tuple[List[GeneMatch], List[GeneMatch]]:
    """
    Split results['regex_findings'] into:
      - temp: records that match any of:
          (a) same gene_id == hgnc_id AND same gene_name == symbol
          (b) same gene_id == hgnc_id
          (c) same gene_name == symbol AND empty gene_id
      - left_to_process: all other records

    Returns (temp, left_to_process).
    """
    # Normalize targets
    """
    temp: records matching any of:
      (a) gene_id == hgnc_id AND gene_name == symbol
      (b) gene_id == hgnc_id
      (c) gene_name == symbol AND gene_id == ""
    left_to_process: all others
    """
    # normalize inputs

    target_symbol = (symbol or "").strip().upper()

    temp: List[GeneMatch] = []
    left: List[GeneMatch] = []

    for gm in results:
        gm_id = (gm.gene_id or "").strip().upper()
        gm_sym = (gm.gene_name or "").strip().upper()

        cond_c = gm_sym == target_symbol and gm_id == ""

        (temp if (cond_c) else left).append(gm)

    return temp, left


def validate_hgnc_ids(hgnc_ids: Set[str]) -> Dict[str, Dict]:
    """Instant O(1) lookups."""
    validated = {}
    for hgnc_id in hgnc_ids:
        # d_num = hgnc_id.replace('HGNC:', '')
        if hgnc_id in HGNC_BY_ID:
            validated[hgnc_id] = HGNC_BY_ID[hgnc_id]
    return validated


def validate_gene_symbols(symbols: Set[str]) -> Dict[str, Dict]:
    """
    Validate a list of potential gene symbols against HGNC_BY_SYMBOL.

    Args:
        symbols: iterable of candidate symbols (e.g., ["TP53", "RRAGD", "foo"])

    Returns:
        Dict mapping the ORIGINAL input symbol -> HGNC row (for approved symbols only).
        Symbols not present in HGNC_BY_SYMBOL are omitted.
    """
    validated: Dict[str, Dict] = {}
    for sym in symbols:
        if not sym:
            continue
        s = sym.strip()
        if not s:
            continue

        # Try exact, then uppercase (common normalization in text)
        row = HGNC_BY_SYMBOL.get(s) or HGNC_BY_SYMBOL.get(s.upper())
        if row:
            validated[sym] = row  # keep original key; value is the HGNC row dict

    return validated


def _process_section_with_tracking(
    section_name: str, section_content: List[str], hgnc_ids_tracker: Set[str]
) -> List[GeneMatch]:
    """Process section and track found HGNC IDs."""

    print(f"\n Processing section: {section_name}")
    section_findings = []

    for element_idx, element_text in enumerate(section_content):
        if not element_text.strip():
            continue

        regex_matches = _extract_with_regex(element_text, section_name, element_idx)

        # Track HGNC IDs for batch validation
        for match in regex_matches:
            if match.gene_id:
                hgnc_ids_tracker.add(match.gene_id)

        section_findings.extend(regex_matches)

    print(f"   Found {len(section_findings)} gene mentions")
    return section_findings


def _extract_with_regex(text: str, section: str, element_idx: int) -> List[GeneMatch]:
    """Extract genes using regex patterns for HGNC format."""
    matches = []

    # Pre-process: Strip parentheses from tokens to catch (HNF1A), (COL4A3), etc.
    # But keep original text for context
    processed_text = re.sub(r"\(([A-Z][A-Z0-9-]{1,11})\)", r"\1", text)

    # Split text into sentences
    sentences = re.split(r"[.!?]+", processed_text)
    original_sentences = re.split(r"[.!?]+", text)

    for sent_idx, (sentence, original_sentence) in enumerate(
        zip(sentences, original_sentences, strict=False)
    ):
        sentence = sentence.strip()
        if not sentence:
            continue

        # PRE-FILTER: Remove stop words from sentence before regex matching
        sentence = _remove_stop_words_from_sentence(sentence)

        #  all findings for this sentence to merge overlapped findings later
        sentence_findings = []

        # Pattern 1: HGNC:1234 format
        hgnc_pattern = r"\bHGNC:(\d+)\b"
        hgnc_matches = re.finditer(hgnc_pattern, sentence)
        for match in hgnc_matches:
            sentence_findings.append(
                GeneMatch(
                    gene_name="",
                    gene_id=f"HGNC:{match.group(1)}",
                    section=section,
                    element_index=element_idx,
                    sentence_index=sent_idx,
                    sentence=original_sentence,
                    method="regex1",
                )
            )

        # Pattern 2: Gene name (HGNC:id) format
        gene_hgnc_pattern = r"\b([A-Z][A-Z0-9-]*)\s*\(\s*HGNC:(\d+)\s*\)"
        gene_hgnc_matches = re.finditer(gene_hgnc_pattern, sentence)
        for match in gene_hgnc_matches:
            sentence_findings.append(
                GeneMatch(
                    gene_name=match.group(1),
                    gene_id=f"HGNC:{match.group(2)}",
                    section=section,
                    element_index=element_idx,
                    sentence_index=sent_idx,
                    sentence=original_sentence,
                    method="regex2",
                )
            )

        # Pattern 4: Gene names with numbers (like COL4A3)
        gene_number_pattern = r"\b([A-Z]{1,11}\d+[A-Z]\d*)\b"
        gene_num_matches = re.finditer(gene_number_pattern, sentence)
        for match in gene_num_matches:
            gene_name = match.group(1)
            if gene_name not in STOP_KEYWORDS:
                sentence_findings.append(
                    GeneMatch(
                        gene_name=gene_name,
                        gene_id="",
                        section=section,
                        element_index=element_idx,
                        sentence_index=sent_idx,
                        sentence=original_sentence,
                        full_element=text,
                        method="regex4",
                    )
                )
        # Pattern 4a: Gene names with numbers (like COL4A3,BRCA1, NPHS2)
        # would not catch BRCA1/2 -need to add later
        gene_number_pattern = r"\b([A-Z]{1,11}\d+(?:[A-Z]\d+)*)\b"
        gene_num_matches = re.finditer(gene_number_pattern, sentence)
        for match in gene_num_matches:
            gene_name = match.group(1)
            if gene_name not in STOP_KEYWORDS:
                sentence_findings.append(
                    GeneMatch(
                        gene_name=gene_name,
                        gene_id="",
                        section=section,
                        element_index=element_idx,
                        sentence_index=sent_idx,
                        sentence=original_sentence,
                        full_element=text,
                        method="regex4a",
                    )
                )

        # Pattern 5: # e.g., RRAGD,
        gene_related_pattern = r"\b([A-Z]{3,11})\b"
        gene_matches = re.finditer(gene_related_pattern, sentence)
        for match in gene_matches:
            gene_name = match.group(1)
            if gene_name != "HGNC" and (gene_name not in STOP_KEYWORDS):
                sentence_findings.append(
                    GeneMatch(
                        gene_name=gene_name,
                        gene_id="",
                        section=section,
                        element_index=element_idx,
                        sentence_index=sent_idx,
                        sentence=original_sentence,
                        full_element=text,
                        method="regex5",
                    )
                )

        # Pattern 3: Gene names with optional -related suffix
        gene_related_pattern = r"\b([A-Z][A-Z0-9]{1,11})(?:-related)\b"
        gene_matches = re.finditer(gene_related_pattern, sentence)
        for match in gene_matches:
            gene_name = match.group(1)
            if gene_name not in STOP_KEYWORDS:
                sentence_findings.append(
                    GeneMatch(
                        gene_name=gene_name,
                        gene_id="",
                        section=section,
                        element_index=element_idx,
                        sentence_index=sent_idx,
                        sentence=original_sentence,
                        full_element=text,
                        method="regex3",
                    )
                )

        # CONSOLIDATION: Collapse overlapping findings
        consolidated_findings = _consolidate_sentence_findings(sentence_findings)

        # Add consolidated findings to main matches
        matches.extend(consolidated_findings)

    return matches


def _consolidate_sentence_findings(findings: List[GeneMatch]) -> List[GeneMatch]:
    """Consolidate overlapping gene findings within a sentence."""
    if not findings:
        return []

    # Group findings by gene_id and gene_name for consolidation
    consolidation_map = {}

    for finding in findings:
        # Create consolidation key
        key = (
            finding.gene_name.upper() if finding.gene_name else "",
            finding.gene_id if finding.gene_id else "",
        )

        if key in consolidation_map:
            # Merge methods
            existing = consolidation_map[key]
            methods = existing.method.split("+") + [finding.method]
            existing.method = "+".join(sorted(set(methods)))
        else:
            consolidation_map[key] = finding

    # Special consolidation: Link gene names with gene IDs
    consolidated = {}

    for (gene_name, gene_id), finding in consolidation_map.items():
        # Look for matches between gene names and gene IDs
        merged = False

        for _existing_key, existing_finding in consolidated.items():
            # Case 1: Same gene ID, merge gene names
            if gene_id and existing_finding.gene_id == gene_id and gene_id != "":
                if gene_name and not existing_finding.gene_name:
                    existing_finding.gene_name = gene_name
                methods = existing_finding.method.split("+") + finding.method.split("+")
                existing_finding.method = "+".join(sorted(set(methods)))
                merged = True
                break

            # Case 2: Same gene name, merge gene IDs
            elif (
                gene_name
                and existing_finding.gene_name.upper() == gene_name.upper()
                and gene_name != ""
            ):
                if gene_id and not existing_finding.gene_id:
                    existing_finding.gene_id = gene_id
                methods = existing_finding.method.split("+") + finding.method.split("+")
                existing_finding.method = "+".join(sorted(set(methods)))
                merged = True
                break

        if not merged:
            consolidated[(gene_name, gene_id)] = finding

    return list(consolidated.values())


def _print_section_results(
    section_name: str, regex_results: List[GeneMatch], total_time: float
) -> None:
    """Print detailed results for a section."""

    print(f"\n📊 SECTION RESULTS: {section_name.upper()}")
    print(f"⏱️  Total processing time: {total_time:.3f} seconds")
    print("-" * 80)

    # Regex results
    print(f"🔤 REGEX FINDINGS ({len(regex_results)} matches):")
    if regex_results:
        for i, match in enumerate(regex_results, 1):
            gene_info = (
                f"{match.gene_name} ({match.gene_id})"
                if match.gene_name and match.gene_id
                else match.gene_id or match.gene_name
            )

            print(f"   {i}. {gene_info} [{match.method}]")
            print(
                f"      Location: Element {match.element_index}, Sentence {match.sentence_index}"
            )
            print(f'      Context: "{match.sentence[:150]}..."')
            print()
    else:
        print("   No regex matches found")

    print("-" * 40)


def _remove_stop_words_from_sentence(sentence: str) -> str:
    """Remove stop word tokens from sentence before gene matching."""

    # Reconstruct sentence with original punctuation/spacing (approximate)
    filtered_sentence = sentence
    for stop_word in STOP_KEYWORDS:
        # Remove stop word as standalone token
        filtered_sentence = re.sub(
            rf"\b{stop_word}\b", "", filtered_sentence, flags=re.IGNORECASE
        )

    return " ".join(filtered_sentence.split())  # Clean up extra spaces


def _print_compressed_findings(all_findings: List[GeneMatch]) -> None:
    """
    DEBUGGING FUNCTION: Print compressed summary of all gene findings before validation.

    Args:
        all_findings: list of all find entities

    Returns:
        None
    """

    if not all_findings:
        print("   No gene findings to display")
        return

    # Group by gene name/ID for compression
    gene_summary = {}

    for finding in all_findings:
        # Create key for grouping
        if finding.gene_name and finding.gene_id:
            key = f"{finding.gene_name} ({finding.gene_id})"
        elif finding.gene_id:
            key = finding.gene_id
        elif finding.gene_name:
            key = finding.gene_name
        else:
            continue

        if key not in gene_summary:
            gene_summary[key] = {"methods": set(), "sections": set(), "count": 0}

        gene_summary[key]["methods"].add(finding.method)
        gene_summary[key]["sections"].add(finding.section)
        gene_summary[key]["count"] += 1

    # Sort by count (most mentioned first)
    sorted_genes = sorted(
        gene_summary.items(), key=lambda x: x[1]["count"], reverse=True
    )

    print(
        f"Found {len(sorted_genes)} unique genes with {len(all_findings)} total mentions:"
    )
    print("-" * 60)

    for i, (gene, info) in enumerate(sorted_genes, 1):
        methods_str = "+".join(sorted(info["methods"]))
        sections_str = ", ".join(sorted(info["sections"]))
        print(f"{i:2d}. {gene}")
        print(
            f"    Methods: [{methods_str}] | Sections: {sections_str} | Mentions: {info['count']}"
        )
        print()
