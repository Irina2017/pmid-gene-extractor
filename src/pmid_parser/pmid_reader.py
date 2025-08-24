"""
pmid_reader: module to get PubMed/PMC content for a single PMID (assessment MVP).

 - Get PubMed metadata via E-utilities from NCBI.
 - Access full-text XML when available; if not - use title + abstract from metadata.
    Out of scope of this MVP: if not full text XML, get free pdf if possible.
 - Return dictionary with metadata and sections content.

NOTE: Assumes paper is in English.
    For MVP, assumes XML has standard sections like Title, Abstract, Introduction
    ( or Background),
    some section about methods, results, discussion, conclusion, references,
    acknowledgements, information about
    authors and supplementary materials.

    Gene mentions in references/supplementary materials are not processed ( out of MVP
    scope).
    Since at least 5 genes are needed, no check in supplementary materials or data in
    this MVP.
    No processing of text information on Figures (out of scope for this MVP due to
    time contraints)
    Tables content would be picked up and processed as text - may require more
    column-aware processing in a future( out of scope of MVP).

Rate limiting unnessary for single PMID processing:
- NCBI limit without API key is 3 req/sec -> only 2 API calls
- If batch processing, consider to use API key later

Performance notes:
- Network latency dominates over XML parsing =>
       2 API calls to get all content, then parsing rather than getting per section
         for parsing.
- One PMID, no batch or asynchronous processing (out of scope) => sequential flow.

Error handling (MVP):
- network timeouts -> exception
- non-existent PMID in database  -> exception [PMID should be validated before -
see helpers.is_valid_pmid()]
- no title or abstract in PubMed -> exception (review XML)
- Network timeouts/invalid PMID -> exception
- Missing title/abstract in PubMed -> exception
- No full PMC access -> warning, metadata-only (title + abstract)
- XML parsing errors -> warning, use of title + abstract
"""

import re
from typing import Dict, List, Optional

import lxml.etree as ET
import requests

from .config import API_TIMEOUT, PMC_API_URL, PUBMED_API_URL


def read_paper_content(pmid: str) -> Dict:
    """
    Get complete paper content for gene extraction.
         Retrieves PubMed metadata (title, abstract, PMC ID, MeSH terms) and PMC
         full-text sections
                  (intro, methods, results, discussion,conclusions) optimized for
                  gene mention analysis.
    Args:
        pmid: Cleaned validated PMID string (e.g., "38790019")

    Returns:
        Dict with 'metadata' (title, abstract, pmc_id,mesh_terms) and 'sections'
        of body
        if free PMC XML available.
        Sections like acknowledgement, references, information about authors and
        supplementary
        materials are skipped for this MVP.
    """
    # Get PMC ID and metadata (first API call)

    pubmed_data = read_pubmed_metadata(pmid)
    print(f"Found article: {pubmed_data['title']}")

    # Get full text if available (second API call)
    if pubmed_data["pmc_id"] is not None:
        print(f"PMC full-text available: {pubmed_data['pmc_id']}")
        pmc_data = read_pmc_fulltext(pubmed_data["pmc_id"])
        sections = extract_sections(pmc_data["raw_xml"])

        # Use existing title/abstract,  should be a match - in metadata

    else:
        print("WARNING: No PMC full-text, using title + abstract only")
        sections = {}

    return {"metadata": pubmed_data, "sections": sections}


def read_pubmed_metadata(pmid: str) -> Dict:
    """
    Get article metadata from PubMed.

    Args:
        pmid: PubMed ID (e.g., "38790019")

    Returns:
        Dict with metadata - pmid, title, abstract, pmc_id, mesh_terms -
        if can retrieve

    Raises:
        requests.RequestException: API call fails
        ValueError: not found PubMed ID
    """

    url = f"{PUBMED_API_URL}?db=pubmed&id={pmid}&retmode=xml"

    try:
        response = requests.get(url, timeout=API_TIMEOUT)
        response.raise_for_status()
    except requests.RequestException as e:
        raise requests.RequestException(f"ERROR: PubMed API failed for PMID {pmid}: {e}") from e

    root = ET.fromstring(response.content)
    if not root.findall("PubmedArticle"):
        raise ValueError(f"ERROR: no such PMID  - {pmid}")

    title = _extract_title(root)
    if not title:
        raise ValueError(f"ERROR: title is not found, inspect XML  - {pmid}")

    abstract = _extract_abstract(root)
    if not abstract:
        raise ValueError(f"ERROR: abstract is not found, inspect XML  - {pmid}")

    pmc_id = _extract_pmc_id(root)
    mesh_terms = _extract_mesh_terms(root)

    return {
        "pmid": pmid,
        "title": title,
        "abstract": abstract,
        "pmc_id": pmc_id,
        "mesh_terms": mesh_terms,  # diseases keywords
    }


def read_pmc_fulltext(pmc_id: str) -> Dict:
    """
    Get full-text XML from PMC.

    Args:
        pmc_id: PMC ID

    Returns:
        Dict with pmc_id and raw xml (as raw_xml)

    Raises:
        requests.RequestException: API call fails

    """
    url = f"{PMC_API_URL}?db=pmc&id={pmc_id}&retmode=xml"

    try:
        response = requests.get(url, timeout=API_TIMEOUT)
        response.raise_for_status()

    except requests.RequestException as e:
        raise requests.RequestException(f"ERROR: PMC API failed for {pmc_id}: {e}") from e

    text = response.content.decode("utf-8", errors="replace")

    # Content-level check: does the XML contain <error>?
    try:
        root = ET.fromstring(text)
        err = root.find(".//error")
        if err is not None:
            # include server’s message for debugging
            raise ValueError(
                f"ERROR: PMC reports error for {pmc_id}: {err.text.strip()
                                                          if err.text else 'Unknown error'}"
            )
    except ET.ParseError as e:
        # Non-XML response (rare), treat as content error
        raise ValueError(f"ERROR: PMC returned non-XML for {pmc_id}") from e

    return {"pmc_id": pmc_id, "raw_xml": text}


def extract_sections(raw_xml: str) -> Dict[str, List[str]]:
    """
    Parse PMC XML into structured sections.

    Args:
        raw_xml: raw xml content

    Returns:
        Dict: {section names (keys) : [text content (values)]}

    """
    try:
        root = ET.fromstring(raw_xml)
    except ET.XMLSyntaxError as e:
        print(f"ERROR: XML parsing failed: {e}, returning empty sections")
        return {}

    # Extract key sections for gene analysis
    sections: Dict[str, List[str]] = {}

    # for abstract and title - will re-use here information
    # from pubmed (should match b/w pubmed and pmc)

    # Introduction: background - > intro (covers introduction )
    #  (only one is available for particular paper)
    sections["intro"] = []
    for section_type in ["background", "intro"]:
        sections["intro"] = _extract_section_text(root, section_type)
        if sections["intro"]:
            break
    print(f"Extracted intro: {len(sections['intro'])}")
    # print(f"Extracted intro: {sections['intro']}")

    # Methods: 'methods', 'materials-methods','materials-and-methods', 'materials_and_methods'
    sections["methods"] = _extract_section_text(root, "methods")
    print(f"Extracted methods: {len(sections['methods'])}")

    # Results: 'results',  'results-discussion','results-and-discussion', 'findings'
    sections["results"] = _extract_section_text(root, "results")
    if not sections["results"]:
        sections["results"] = _extract_section_text(root, "findings")
    print(f"Extracted results: {len(sections['results'])}")

    # for p in sections['results']:
    #    print("****")
    #    print(len(p))
    #    print(p)
    #    print("****")

    # Only discussion here - if merged with results would be picked up before
    sections["discussion"] = _extract_section_text(root, "discussion")
    print(f"Extracted discussion: {len(sections['discussion'])}")

    sections["conclusions"] = _extract_section_text(root, "conclusion")  # conclusions
    print(f"Extracted conclusions: {len(sections['conclusions'])}")

    # Mentions of genes in references are not counted as mentions in the current PMID (MVP).
    # skipping 'references', 'acknowledgments','author-information', 'supplementary-material'

    # Remove empty sections
    sections = {k: v for k, v in sections.items() if v and any(p.strip() for p in v)}

    return sections


# Private functions for XML parsing
def _extract_title(root: ET.Element) -> str:
    """
    Extract article title from PubMed XML

    Args:
        root: XML root element


    Returns:
        string with title or empty if not found
    """

    title_elem = root.find(".//ArticleTitle")

    if title_elem is not None:
        # Get all text content, including from nested elements
        title_text = "".join(title_elem.itertext()).strip()

        if title_text:
            return title_text

    return ""


def _extract_abstract(root: ET.Element) -> str:
    """Extract abstract text from PubMed XML.

    Args:
        root: XML root element


    Returns:
        Combined abstract text
    """

    elements = root.findall(".//AbstractText")
    if elements:
        # need to account for <i> and <b> etc
        texts = ["".join(elem.itertext()).strip() for elem in elements]
        result = " ".join(t for t in texts if t)
        if result:
            return result
    return ""


def _extract_pmc_id(root: ET.Element) -> Optional[str]:
    """Extract PMC ID if available."""
    for article_id in root.findall(".//ArticleId"):
        if article_id.get("IdType") == "pmc":
            return article_id.text
    return None


def _extract_mesh_terms(root: ET.Element) -> List[str]:
    """Extract MeSH disease terms."""
    mesh_terms = []
    for descriptor in root.findall(".//DescriptorName"):
        if descriptor.text:
            mesh_terms.append(descriptor.text)
    return mesh_terms


def _extract_section_text(root: ET.Element, section_type: str) -> List[str]:
    """Extract text from specific PMC body section

    Args:
        root: PMC XML root element
        section_type: section to extract

    Returns:
        List of paragraph texts from the section, or empty list if not found

    """

    # try sec-type attribute - if lucky for faster processing
    paragraphs = root.findall(f'.//sec[@sec-type="{section_type}"]//p')
    result = _combine_paragraphs(paragraphs)
    if result:
        return result

    # data-title attribute is not consistenly present
    # go with <title> - required part of <sec>
    # try <title> text content (case-insensitive contains)

    xpath_query = (
        f".//sec[title[contains("
        f'translate(text(), "ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz"), '
        f'"{section_type.lower()}")]]//p'
    )
    paragraphs = root.xpath(xpath_query)

    # Filter out nested paragraphs for results section -caused by Results in Figures or Tables titles
    if section_type == "results":
        paragraphs = _filter_nested_paragraphs(paragraphs)
        print(f"After filtering nested: {len(paragraphs)} paragraphs")

    # if section_type=="results":
    #    for i, p in enumerate(paragraphs):
    #        print(f"--- Paragraph {i+1} ---")
    #        print(f"Tag: {p.tag}")
    #        print(f"Has children: {len(list(p)) > 0}")
    #        print(f"Full content: {' '.join(p.itertext()).strip()}")
    #        print()

    result = _combine_paragraphs(paragraphs)
    if result:
        return result

    return []


def _combine_paragraphs(paragraphs: List[ET.Element]) -> List[str]:
    """
    Combine text from list of paragraph elements.

    Args:
        paragraphs: List of paragraph XML elements

    Returns:
        List of paragraph texts (non-empty only)
    """
    if paragraphs:
        texts = []
        for p in paragraphs:
            full_text = " ".join(p.itertext()).strip()
            # if full_text:
            #    texts.append(full_text)
            if full_text:
                # Remove citations
                clean_text = re.sub(r"\[\s*\d+(?:\s*[-–,]\s*\d+)*\s*\]", "", full_text)

                # Clean up extra whitespace from removed citations
                clean_text = " ".join(clean_text.split())

                if clean_text:
                    texts.append(clean_text)
        return texts
    return []


def _filter_nested_paragraphs(paragraphs: List[ET.Element]) -> List[ET.Element]:
    """Remove paragraphs that are nested inside other paragraphs in the same list.

    Args:
        paragraphs: List of paragraph elements

    Returns:
        Filtered list with nested paragraphs removed
    """
    if not paragraphs:
        return []

    filtered_paragraphs = []

    for p in paragraphs:
        is_nested = False

        # Check if this paragraph is a child of any other paragraph in the list
        for other_p in paragraphs:
            if p != other_p:
                if _is_descendant(p, other_p):
                    is_nested = True
                    break

        # Only keep paragraphs that are not nested inside others
        if not is_nested:
            filtered_paragraphs.append(p)

    return filtered_paragraphs


def _is_descendant(child: ET.Element, potential_parent: ET.Element) -> bool:
    """
    Check if child element is a descendant of potential_parent.

    Args:
          child: one element in XML
          potential_parent: another element in XML

    Returns:
            True if nested inside, false - otherwise


    """
    current = child.getparent()
    while current is not None:
        if current == potential_parent:
            return True
        current = current.getparent()
    return False
