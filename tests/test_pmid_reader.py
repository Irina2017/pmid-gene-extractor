"""
pytest tests for pmid_reader module.

Not exhaustive due to time constraints.
"""

import time
from pathlib import Path

import lxml.etree as ET
import pytest
import requests

from pmid_parser.config import RATE_LIMIT_DELAY_SECONDS

# Assuming your module is named pmid_reader
from pmid_parser.pmid_reader import (
    _combine_paragraphs,
    _extract_section_text,
    _filter_nested_paragraphs,
    _is_descendant,
    read_paper_content,
    read_pmc_fulltext,
    read_pubmed_metadata,
)


# to avoid API problems with PubMed or PMC
@pytest.fixture(autouse=True)
def rate_guard():
    time.sleep(RATE_LIMIT_DELAY_SECONDS)


REAL_PUBMED_TEST_DATA = [
    {
        "pmid": "38790019",
        "title": "Diagnostic yield of exome and genome sequencing after non-diagnostic "
        "multi-gene panels in patients with single-system diseases.",
        "abstract_parts": [
            "Though next-generation sequencing (NGS) tests like exome sequencing (ES), genome "
            "sequencing (GS), and panels derived from exome and genome data (EGBP) are "
            "effective for rare diseases, the ideal diagnostic approach is debated. Limited "
            "research has explored reanalyzing raw ES and GS data post-negative EGBP results for diagnostics.",
            "We analyzed complete ES/GS raw sequencing data from Mayo Clinic's Program for "
            "Rare and Undiagnosed Diseases (PRaUD) patients to assess whether supplementary "
            "findings could augment diagnostic yield. ES data from 80 patients (59 adults) "
            "and GS data from 20 patients (10 adults), averaging 43 years in age, were analyzed."
            " Most patients had renal (n=44) and auto-inflammatory (n=29) phenotypes. "
            "Ninety-six cases had negative findings and in four cases additional genetic variants "
            "were found, including a variant related to a recently described disease (RRAGD-related"
            " hypomagnesemia), a variant missed due to discordant inheritance pattern (COL4A3), "
            "a variant with high allelic frequency (NPHS2) in the general population, and a variant "
            "associated with an initially untargeted phenotype (HNF1A).",
            "ES and GS show diagnostic yields comparable to EGBP for single-system diseases. "
            "However, EGBP's limitations in detecting new disease-associated genes underscore the "
            "necessity for periodic updates.",
        ],
        "pmc_id": "PMC11127317",
        "mesh_terms": [
            "Humans",
            "Adult",
            "Female",
            "Male",
            "Middle Aged",
            "High-Throughput Nucleotide Sequencing",
            "Exome Sequencing",
            "Exome",
            "Young Adult",
            "Rare Diseases",
            "Aged",
            "Adolescent",
            "Whole Genome Sequencing",
        ],
    },
    {
        "pmid": "38674358",
        "title": "Expansion of the Genotypic and Phenotypic Spectrum of ASH1L-Related Syndromic "
        "Neurodevelopmental Disorder.",
        "abstract_parts": [
            "Pathogenic ASH1L variants have been reported in probands with broad phenotypic "
            "presentations, including intellectual disability, autism spectrum disorder, attention "
            "deficit hyperactivity disorder, seizures, congenital anomalies, and other skeletal, muscular, "
            "and sleep differences. Here, we review previously published individuals with pathogenic ASH1L "
            "variants and report three further probands with novel ASH1L variants and previously unreported "
            "phenotypic features, including mixed receptive language disorder and gait disturbances. "
            "These novel data from the Brain Gene Registry, an accessible repository of clinically "
            "derived genotypic and phenotypic data, have allowed for the expansion of the phenotypic"
            " and genotypic spectrum of this condition."
        ],
        "pmc_id": "PMC11049257",
        "mesh_terms": [
            "Humans",
            "Neurodevelopmental Disorders",
            "Phenotype",
            "Male",
            "Histone-Lysine N-Methyltransferase",
            "Female",
            "Child",
            "Genotype",
            "DNA-Binding Proteins",
            "Intellectual Disability",
            "Transcription Factors",
            "Child, Preschool",
            "Autism Spectrum Disorder",
            "Mutation",
            "Adolescent",
        ],
    },
    {
        "pmid": "39442041",
        "title": "Neurodevelopmental Disorder Caused by Deletion of CHASERR,"
        " a lncRNA Gene.",
        "abstract_parts": [
            "CHASERR encodes a human long noncoding RNA (lncRNA) adjacent to CHD2, "
            "a coding gene in which de novo loss-of-function variants cause developmental "
            "and epileptic encephalopathy. Here, we report our findings in three "
            "unrelated children with a syndromic, early-onset neurodevelopmental disorder, "
            "each of whom had a de novo deletion in the CHASERR locus. The children had severe "
            "encephalopathy, shared facial dysmorphisms, cortical atrophy, and cerebral hypomyelination"
            " - a phenotype that is distinct from the phenotypes of patients with CHD2 haploinsufficiency."
            " We found that the CHASERR deletion results in increased CHD2 protein abundance in "
            "patient-derived cell lines and increased expression of the CHD2 transcript in cis. "
            "These findings indicate that CHD2 has bidirectional dosage sensitivity in human disease, "
            "and we recommend that other lncRNA-encoding genes be evaluated, particularly those upstream of "
            "genes associated with mendelian disorders. (Funded by the National Human Genome Research"
            " Institute and others.)."
        ],
        "pmc_id": "PMC11826417",
        "mesh_terms": [
            "Child, Preschool",
            "Female",
            "Humans",
            "Infant",
            "Male",
            "Brain",
            "DNA-Binding Proteins",
            "Gene Deletion",
            "Haploinsufficiency",
            "Neurodevelopmental Disorders",
            "Phenotype",
            "RNA, Long Noncoding",
            "Sequence Deletion",
        ],
    },
    {
        "pmid": "16804717",
        "title": "Laryngeal leiomyosarcoma: a case report and review of the literature.",
        "abstract_parts": [
            "Laryngeal leiomyosarcoma (LLM) is a rare malignancy originating "
            "from the smooth muscles of blood vessels or from aberrant "
            "undifferentiated mesenchymal tissue. Histological diagnosis"
            " may be particularly difficult and correct diagnosis "
            "is based on immunohistochemical investigations and electron microscopy. "
            "A case report of a LLM in a 74-year-old man is presented. Direct laryngoscopy "
            "revealed a large glottic lesion causing airway compromise and an emergency "
            "tracheotomy was performed. Subsequent total laryngectomy confirmed the"
            " diagnosis of leiomyosarcoma. Lung metastases developed 8 months following "
            "treatment, despite the absence of local or regional recurrence, and the "
            "patient died 3 months later. A review of the English and French "
            "literature revealed 30 previous cases of LLM. Clinical presentation, "
            "histological diagnosis, and management of this rare malignancy "
            "are analyzed aiming to improve our knowledge regarding the best treatment modality."
        ],
        "pmc_id": None,  # No free access
        "mesh_terms": [
            "Aged",
            "Diagnosis, Differential",
            "Fatal Outcome",
            "Humans",
            "Laryngeal Neoplasms",
            "Laryngectomy",
            "Laryngoscopy",
            "Leiomyosarcoma",
            "Lung Neoplasms",
            "Male",
        ],
    },
]


# due to time constraint, test helpers functions with read_pubmed_metadata
@pytest.mark.parametrize("expected", REAL_PUBMED_TEST_DATA)
def test_read_pubmed_metadata_real_pmids(expected):
    """Integration test for read_pubmed_metadata with REAL PMID data ."""

    # Make actual API call to PubMed
    result = read_pubmed_metadata(expected["pmid"])

    print(f"\n=== PMID {expected['pmid']} ===")

    expected_abstract = " ".join(expected["abstract_parts"])
    # Assert results match expected data
    assert result["pmc_id"] == expected["pmc_id"]
    assert result["title"] == expected["title"]
    assert result["abstract"] == expected_abstract
    assert result["mesh_terms"] == expected["mesh_terms"]


@pytest.mark.parametrize(
    "expected",
    [
        {
            "metadata": {
                "pmid": "16804717",
                "title": "Laryngeal leiomyosarcoma: a case report and review of the "
                "literature.",
                "abstract_parts": [
                    "Laryngeal leiomyosarcoma (LLM) is a rare malignancy originating "
                    "from the smooth muscles of blood vessels or from aberrant undifferentiated"
                    " mesenchymal tissue. "
                    "Histological diagnosis may be particularly difficult and correct diagnosis "
                    "is based on immunohistochemical investigations and electron microscopy."
                    " A case report of a LLM in a 74-year-old man is presented. Direct "
                    "laryngoscopy revealed a large glottic lesion causing airway compromise"
                    " and an emergency tracheotomy was performed. Subsequent total "
                    "laryngectomy confirmed the diagnosis of leiomyosarcoma."
                    " Lung metastases developed 8 months following treatment, "
                    "despite the absence of local or regional recurrence, and"
                    " the patient died 3 months later. A review of the English "
                    "and French literature revealed 30 previous cases of LLM."
                    " Clinical presentation, histological diagnosis, and "
                    "management of this rare malignancy are analyzed aiming"
                    " to improve our knowledge regarding the best treatment modality."
                ],
                "pmc_id": None,  # No free access
                "mesh_terms": [
                    "Aged",
                    "Diagnosis, Differential",
                    "Fatal Outcome",
                    "Humans",
                    "Laryngeal Neoplasms",
                    "Laryngectomy",
                    "Laryngoscopy",
                    "Leiomyosarcoma",
                    "Lung Neoplasms",
                    "Male",
                ],
            },
            "sections": {},
        }
    ],
)
def test_read_paper_content_no_free_access(expected):

    result = read_paper_content(expected["metadata"]["pmid"])

    print(f"\n=== PMID {expected['metadata']["pmid"]} ===")

    expected_abstract = " ".join(expected["metadata"]["abstract_parts"])
    # Assert results match expected data
    assert result["metadata"]["pmc_id"] == expected["metadata"]["pmc_id"]
    assert result["metadata"]["title"] == expected["metadata"]["title"]
    assert result["metadata"]["abstract"] == expected_abstract
    assert result["metadata"]["mesh_terms"] == expected["metadata"]["mesh_terms"]
    assert result["sections"] == {}


def test_read_pubmed_metadata_fake_pmid():
    # was not found at Aug 23, 2025
    with pytest.raises(ValueError) as error_info:
        read_pubmed_metadata("98790018")
    assert "ERROR: no such PMID  - 98790018" in str(error_info)


def test_read_pubmed_metadata_empty_pmid():
    with pytest.raises(requests.exceptions.RequestException) as error_info:
        read_pubmed_metadata("pmid")  # not possible to get here through CLI
    assert (
        "ERROR: PubMed API failed for PMID pmid: 400 Client Error: Bad Request for url:"
        in str(error_info)
    )


XML = """
<root>
  <!-- will produce "Alpha BRCA1 signal." -->
  <p id="p1">Alpha<italic>BRCA1</italic>signal.</p>
  <!-- empty / whitespace-only -->
  <p id="p2"></p>
  <p id="p5">   </p>
  <p id="p3">Outer<b>bold</b>and<i>italic</i>text</p>
  <p id="p4"><xref>Figure 1</xref>:Details.</p>
</root>
"""


def test_combine_paragraphs_filters_empty_and_preserves_order():
    root = ET.fromstring(XML)
    p1 = root.xpath(".//p[@id='p1']")[0]
    p2 = root.xpath(".//p[@id='p2']")[0]
    p3 = root.xpath(".//p[@id='p3']")[0]
    p4 = root.xpath(".//p[@id='p4']")[0]
    p5 = root.xpath(".//p[@id='p5']")[0]

    # Include empty/whitespace-only paragraphs; function should skip them
    out = _combine_paragraphs([p1, p2, p3, p4, p5])

    assert out == [
        "Alpha BRCA1 signal.",
        "Outer bold and italic text",
        "Figure 1 :Details.",
    ]


def test_combine_paragraphs_empty_input_returns_empty_list():
    assert _combine_paragraphs([]) == []


def test_combine_paragraphs_removes_citations():
    """Test that citation patterns are removed from paragraph text."""
    xml_content = """
    <root>
        <p>This text has citation [13] in middle.</p>
        <p>Multiple citations [ 18 , 19 ] with spaces.</p>
        <p>Range citation [5-8] and list [1,2,3] citations.</p>
        <p>Text with [42] and normal text after.</p>
    </root>
    """
    root = ET.fromstring(xml_content)
    paragraphs = root.findall(".//p")

    result = _combine_paragraphs(paragraphs)

    assert result == [
        "This text has citation in middle.",
        "Multiple citations with spaces.",
        "Range citation and list citations.",
        "Text with and normal text after.",
    ]


def test_combine_paragraphs_preserves_figure_table_ref():
    """Test that Figure/Table references are ok after removal of citations."""
    xml_content = """
    <root>
        <p>See <a href="#fig1">Figure 1</a> for details [13].</p>
        <p>Results in <a href="#table2">Table 2</a> show significance [ 18 , 19 ].</p>
        <p>Refer to <a href="#supp">Supplementary Table S1</a> data [5-8].</p>
    </root>
    """
    root = ET.fromstring(xml_content)
    paragraphs = root.findall(".//p")

    result = _combine_paragraphs(paragraphs)

    assert result == [
        # the extra space is fine for the purpose of this project (another re to clean it later)
        "See Figure 1 for details .",
        "Results in Table 2 show significance .",
        "Refer to Supplementary Table S1 data .",
    ]


def test_combine_paragraphs_empty_after_citation_removal():
    """Test paragraph that becomes empty after citation removal."""
    xml_content = """
    <root>
        <p>[13, 14, 15]</p>
        <p>   [ 42 ]   </p>
        <p>Real content here.</p>
        <p></p>
    </root>
    """
    root = ET.fromstring(xml_content)
    paragraphs = root.findall(".//p")

    result = _combine_paragraphs(paragraphs)

    # Only the paragraph with real content should remain
    assert result == ["Real content here."]


XML2 = """
<root>
  <sec id="A">
    <p id="outer">
      Outer
      <p id="inner">Inner</p>
    </p>
    <p id="sibA">Sibling A</p>
  </sec>
  <sec id="B">
    <p id="other">Other</p>
  </sec>
</root>
"""


def test_filter_nested_paragraphs_removes_only_nested_ps():
    root = ET.fromstring(XML2)
    outer = root.xpath(".//p[@id='outer']")[0]
    inner = root.xpath(".//p[@id='inner']")[0]
    sibA = root.xpath(".//p[@id='sibA']")[0]
    other = root.xpath(".//p[@id='other']")[0]

    paragraphs = [outer, inner, sibA, other]
    filtered = _filter_nested_paragraphs(paragraphs)

    # inner should be removed; order preserved for remaining
    assert [p.get("id") for p in filtered] == ["outer", "sibA", "other"]


def test_filter_nested_paragraphs_no_parent():
    """If the parent isn't in the list, a nested <p> should NOT be removed."""
    root = ET.fromstring(XML2)
    inner = root.xpath(".//p[@id='inner']")[0]
    other = root.xpath(".//p[@id='other']")[0]

    paragraphs = [inner, other]
    filtered = _filter_nested_paragraphs(paragraphs)

    assert [p.get("id") for p in filtered] == ["inner", "other"]


def test_is_descendant_p_inside_p():

    root = ET.fromstring(XML2)

    outer = root.xpath(".//p[@id='outer']")[0]
    inner = root.xpath(".//p[@id='inner']")[0]
    grandparent = root.xpath(".//sec[@id='A']")[0]
    sibA = root.xpath(".//p[@id='sibA']")[0]
    other = root.xpath(".//p[@id='other']")[0]

    # child → parent
    assert _is_descendant(inner, outer) is True
    # child → grandparent (sec A)
    assert _is_descendant(inner, grandparent) is True
    # siblings are not nested
    assert _is_descendant(sibA, other) is False


def create_test_xml():
    """Create realistic PMC XML with various section structures and nested elements."""
    xml_content = """
    <root>
          <sec sec-type="intro">
          <title>Introduction</title>
            <p>Introduction paragraph 1.</p>
            <p>Introduction paragraph 2.</p>
            <p>Introduction paragraph 3 - Table 1.</p>
         </sec>
       <sec>
            <title>Materials and Methods</title>
            <p>Methods paragraph 1 with gene <a href="#ref1">BRCA1</a> mentioned.</p>
            <p>Methods paragraph 2 has nested content: Sub-item 1 Sub-item 2 and
            external <a href="http://example.com">database link</a>.</p>
            <p>Methods paragraph 3 is clean text with <a href="#fig1">Figure 1</a>
            reference.</p>
        </sec>
        <sec>
          <title>Results</title>
           <p>Main results paragraph 1 shows significant findings.</p>
           <p>Results paragraph 2 references <a href="#fig1">Figure 1A </a> and
           <a href="#table1">Table 2</a>.</p>
         <p>Results paragraph 3 after nested section mentions <a href="#fig2">
         Figure 2B</a>.</p>
         <p>Final results paragraph discusses <a href="#table2">Supplementary
         Table S1</a> data.</p>
        </sec>
         <sec>
          <title>DISCUSSION</title>
           <p>Main points.</p>
            <p></p>
         </sec>
</root>
    """
    return xml_content


def test_extract_by_sec_type_attribute_intro():

    xml_content = create_test_xml()
    root = ET.fromstring(xml_content)

    # Test introduction section
    result = _extract_section_text(root, "intro")

    assert len(result) == 3
    print(result)

    assert "Introduction paragraph 1." in result[0]
    assert "Introduction paragraph 2." in result[1]
    assert "Table 1" in result[2]


def test_extract_methods_section():
    """Test extraction of methods section with hyperlinks and references."""
    xml_content = create_test_xml()
    root = ET.fromstring(xml_content)

    # Test methods section
    result = _extract_section_text(root, "methods")

    assert len(result) == 3

    assert result == [
        "Methods paragraph 1 with gene BRCA1 mentioned.",
        "Methods paragraph 2 has nested content: Sub-item 1 Sub-item 2 and "
        "external database link .",
        "Methods paragraph 3 is clean text with Figure 1 reference.",
    ]


def test_extract_results_section():
    """Test extraction of results section"""
    xml_content = create_test_xml()
    root = ET.fromstring(xml_content)

    # Test methods section
    result = _extract_section_text(root, "results")

    assert len(result) == 4

    assert result == [
        "Main results paragraph 1 shows significant findings.",
        "Results paragraph 2 references Figure 1A and Table 2 .",
        "Results paragraph 3 after nested section mentions Figure 2B .",
        "Final results paragraph discusses Supplementary Table S1 data.",
    ]


def test_extract_discussion_section():
    """Test extraction of discussion section"""
    xml_content = create_test_xml()
    root = ET.fromstring(xml_content)

    # Test methods section
    result = _extract_section_text(root, "discussion")

    assert len(result) == 1
    assert result == ["Main points."]


def test_extract_section_not_in_xml():
    """Test extraction of discussion section"""
    xml_content = create_test_xml()
    root = ET.fromstring(xml_content)

    # Test methods section
    result = _extract_section_text(root, "something")
    assert len(result) == 0
    assert result == []


def test_read_pmc_fulltext():
    expected = (Path(__file__).parent / "data" / "PMC11127315.xml").read_text(
        encoding="utf-8"
    )
    result = read_pmc_fulltext("PMC11127315")
    assert result["raw_xml"].strip() == expected.strip()


def test_read_pmc_fulltext_fake_pmid():
    # was not found at Aug 23, 2025
    with pytest.raises(ValueError) as error_info:
        read_pmc_fulltext("PMC91127315")
    assert (
        "ERROR: PMC reports error for PMC91127315: "
        "The following PMCID is not available: 91127315" in str(error_info)
    )
