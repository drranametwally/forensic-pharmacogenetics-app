import io
import re
import csv
import hashlib
from datetime import datetime
from typing import Dict, List, Tuple, Optional

import pandas as pd
import streamlit as st

# Optional PDF dependency
try:
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
    )
    REPORTLAB_OK = True
except Exception:
    REPORTLAB_OK = False


# ============================================================
# R-FINGERPRINT DNA SUITE
# FORENSIC PHARMACOGENETICS & TOXICOLOGY DECISION SUPPORT
# Version 3.0 — Validation-Ready Research Prototype
# ============================================================
#
# IMPORTANT:
# This application is a research/decision-support prototype.
# It is NOT independently validated as an IVD, diagnostic device,
# or autonomous clinical decision-maker. Local analytical and
# clinical validation, laboratory SOPs, privacy controls, ethics
# approval, and professional review are required before real-world
# clinical deployment.
#
# Core pipeline:
# Input -> QC -> Variant parsing -> Genotype evidence ->
# Phenotype interpretation -> Drug/substance-specific evidence ->
# Toxicology risk indicator -> Explainable report
#
# The engine deliberately avoids assigning "high risk" merely because
# a file is missing. Missing/unsupported data produce an explicit
# insufficient-data or review-required status.
# ============================================================

st.set_page_config(
    page_title="R-Fingerprint DNA Suite",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

APP_VERSION = "3.0"
BUILD_DATE = datetime.now().strftime("%Y-%m-%d")

GENES = [
    "CYP2D6", "CYP2C19", "CYP3A4", "OPRM1", "ABCB1",
    "CYP2B6", "CYP2C9", "HTR2A", "ANKK1", "COMT", "SLCO1B1", "HLA-B"
]

# Variant-level mappings intentionally remain conservative.
# A single SNV is NOT treated as a complete star-allele call unless
# the user supplies a supported diplotype or a validated upstream
# pharmacogenomics caller.
VARIANT_MAP = {
    "rs3892097": {"gene": "CYP2D6", "label": "CYP2D6*4-defining variant", "allele": "*4"},
    "rs35742686": {"gene": "CYP2D6", "label": "CYP2D6*3-defining variant", "allele": "*3"},
    "rs1065852": {"gene": "CYP2D6", "label": "CYP2D6*10-associated variant", "allele": "*10"},
    "rs5030865": {"gene": "CYP2D6", "label": "CYP2D6*14-associated variant", "allele": "*14"},
    "rs1799971": {"gene": "OPRM1", "label": "OPRM1 A118G", "allele": "A118G"},
    "rs4244285": {"gene": "CYP2C19", "label": "CYP2C19*2", "allele": "*2"},
    "rs12248560": {"gene": "CYP2C19", "label": "CYP2C19*17-associated variant", "allele": "*17"},
    "rs28399504": {"gene": "CYP2C19", "label": "CYP2C19*2-associated variant", "allele": "*2"},
    "rs35599367": {"gene": "CYP3A4", "label": "CYP3A4*22-associated variant", "allele": "*22"},
    "rs4149056": {"gene": "SLCO1B1", "label": "SLCO1B1 c.521T>C", "allele": "*5/*15-associated"},
    "rs1045642": {"gene": "ABCB1", "label": "ABCB1 3435C>T", "allele": "3435C>T"},
    "rs4149056": {"gene": "SLCO1B1", "label": "SLCO1B1 c.521T>C", "allele": "*5/*15-associated"},
    "rs2032582": {"gene": "ABCB1", "label": "ABCB1 rs2032582 (G2677T/A)", "allele": "transporter variant", "category": "Drug transport", "evidence_level": "Research / drug-specific", "clinical_actionability": False},
    "rs1799853": {"gene": "CYP2C9", "label": "CYP2C9*2-defining variant", "allele": "*2", "category": "Drug metabolism", "evidence_level": "Established for relevant drugs", "clinical_actionability": True},
    "rs1057910": {"gene": "CYP2C9", "label": "CYP2C9*3-defining variant", "allele": "*3", "category": "Drug metabolism", "evidence_level": "Established for relevant drugs", "clinical_actionability": True},
    "rs1800497": {"gene": "ANKK1", "label": "ANKK1 rs1800497 (Taq1A; historically discussed as DRD2-associated)", "allele": "Taq1A", "category": "Substance-use susceptibility", "evidence_level": "Research / inconsistent", "clinical_actionability": False},
    "rs6311": {"gene": "HTR2A", "label": "HTR2A rs6311", "allele": "-1438G>A", "category": "Psychotropic response", "evidence_level": "Research / drug-specific", "clinical_actionability": False},
    "rs6313": {"gene": "HTR2A", "label": "HTR2A rs6313 (T102C)", "allele": "T102C", "category": "Psychotropic response", "evidence_level": "Research / drug-specific", "clinical_actionability": False},
    "rs4680": {"gene": "COMT", "label": "COMT Val158Met", "allele": "Val158Met", "category": "Pain / opioid response", "evidence_level": "Mixed / limited", "clinical_actionability": False},
    "rs2395029": {"gene": "HLA-B", "label": "HLA-linked marker; requires validated HLA interpretation", "allele": "marker"},
}

# Conservative, drug-specific rules. These are intentionally not a universal
# "toxicity score". They express documented pharmacogenetic concerns.
DRUG_RULES = {
    "Codeine": [
        {
            "gene": "CYP2D6",
            "phenotypes": ["Ultrarapid metabolizer"],
            "indicator": "Increased toxicity concern",
            "severity": "High",
            "mechanism": "More rapid conversion of codeine to morphine can increase morphine exposure.",
            "evidence": "CPIC CYP2D6-opioid guidance",
            "source": "https://cpicpgx.org/gene/CYP2D6/",
        },
        {
            "gene": "CYP2D6",
            "phenotypes": ["Poor metabolizer"],
            "indicator": "Reduced analgesic effect concern",
            "severity": "Moderate",
            "mechanism": "Reduced conversion of codeine to morphine may reduce analgesic effect; this is not equivalent to proven overdose risk.",
            "evidence": "CPIC CYP2D6-opioid guidance",
            "source": "https://cpicpgx.org/gene/CYP2D6/",
        },
        {
            "gene": "CYP2D6",
            "phenotypes": ["Intermediate metabolizer"],
            "indicator": "Reduced/variable response concern",
            "severity": "Moderate",
            "mechanism": "Reduced CYP2D6 activity can reduce formation of active metabolite.",
            "evidence": "CPIC CYP2D6-opioid guidance",
            "source": "https://cpicpgx.org/gene/CYP2D6/",
        },
    ],
    "Tramadol": [
        {
            "gene": "CYP2D6",
            "phenotypes": ["Ultrarapid metabolizer"],
            "indicator": "Increased toxicity concern",
            "severity": "High",
            "mechanism": "Increased formation of active metabolite may increase toxicity risk.",
            "evidence": "CPIC CYP2D6-opioid guidance",
            "source": "https://cpicpgx.org/gene/CYP2D6/",
        },
        {
            "gene": "CYP2D6",
            "phenotypes": ["Poor metabolizer"],
            "indicator": "Reduced analgesic effect concern",
            "severity": "Moderate",
            "mechanism": "Reduced formation of active metabolite can reduce analgesic effect.",
            "evidence": "CPIC CYP2D6-opioid guidance",
            "source": "https://cpicpgx.org/gene/CYP2D6/",
        },
    ],
    "Clopidogrel": [
        {
            "gene": "CYP2C19",
            "phenotypes": ["Intermediate metabolizer", "Poor metabolizer"],
            "indicator": "Reduced activation / reduced platelet inhibition concern",
            "severity": "High",
            "mechanism": "CYP2C19 contributes to clopidogrel bioactivation; reduced function can reduce active metabolite formation.",
            "evidence": "CPIC CYP2C19-clopidogrel guideline",
            "source": "https://files.cpicpgx.org/data/guideline/publication/clopidogrel/2022/35034351.pdf",
        }
    ],
    "Amitriptyline": [
        {
            "gene": "CYP2D6",
            "phenotypes": ["Ultrarapid metabolizer", "Poor metabolizer"],
            "indicator": "Exposure/response concern",
            "severity": "Moderate",
            "mechanism": "CYP2D6 phenotype may alter systemic concentrations; interpretation should follow the applicable guideline/label.",
            "evidence": "FDA pharmacogenetic association table",
            "source": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
        },
        {
            "gene": "CYP2C19",
            "phenotypes": ["Poor metabolizer", "Intermediate metabolizer"],
            "indicator": "Exposure/response concern",
            "severity": "Moderate",
            "mechanism": "CYP2C19 phenotype can alter exposure; interpretation is drug-specific.",
            "evidence": "FDA pharmacogenetic association table",
            "source": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
        },
    ],
    "Omeprazole": [
        {
            "gene": "CYP2C19",
            "phenotypes": ["Poor metabolizer", "Intermediate metabolizer"],
            "indicator": "Increased exposure concern",
            "severity": "Moderate",
            "mechanism": "Reduced CYP2C19 activity can increase systemic exposure.",
            "evidence": "FDA pharmacogenetic association table",
            "source": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
        }
    ],
    "Methadone": [
        {
            "gene": "CYP2B6",
            "phenotypes": ["Reduced function", "Poor metabolizer"],
            "indicator": "Exposure/response concern",
            "severity": "Moderate",
            "mechanism": "CYP2B6 pharmacogenetic variation can influence methadone disposition; interpretation requires drug-specific validated evidence.",
            "evidence": "CPIC methadone/CYP2B6 guideline resources",
            "source": "https://cpicpgx.org/gene/CYP2B6/",
        }
    ],
}

SEVERITY_RANK = {"Low": 1, "Moderate": 2, "High": 3, "Review": 4}

# -----------------------------
# Styling
# -----------------------------
st.markdown("""
<style>
.main-title {font-size: 34px; font-weight: 800; margin-bottom: 4px;}
.subtitle {font-size: 16px; color: #666; margin-bottom: 20px;}
.section-title {font-size: 22px; font-weight: 750; margin-top: 18px;}
.small-muted {font-size: 12px; color: #666;}
.status-box {padding: 12px; border-radius: 10px; border: 1px solid #ddd;}
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-title">🧬 R-Fingerprint DNA Suite</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="subtitle">Forensic Pharmacogenetics & Toxicology Decision-Support — Validation-Ready Research Prototype</div>',
    unsafe_allow_html=True
)

# -----------------------------
# Utility functions
# -----------------------------
def normalize_gene(value: str) -> str:
    value = str(value).strip().upper().replace(" ", "")
    aliases = {
        "CYP2D6": "CYP2D6",
        "CYP2C19": "CYP2C19",
        "CYP3A4": "CYP3A4",
        "OPRM1": "OPRM1",
        "ABCB1": "ABCB1",
        "CYP2B6": "CYP2B6",
        "CYP2C9": "CYP2C9",
        "HTR2A": "HTR2A",
        "ANKK1": "ANKK1",
        "DRD2": "ANKK1",
        "COMT": "COMT",
        "SLCO1B1": "SLCO1B1",
        "HLA-B": "HLA-B",
        "HLAB": "HLA-B",
    }
    return aliases.get(value, value)


def normalize_rs(value: str) -> str:
    value = str(value).strip().lower()
    if value.startswith("rs"):
        return value
    m = re.search(r"(rs\d+)", value)
    return m.group(1) if m else value


def genotype_from_gt(gt: str) -> str:
    if gt is None:
        return ""
    gt = str(gt).replace("|", "/")
    return gt.split(":")[0]


def allele_pair_from_gt(gt: str, ref: str, alt: str) -> str:
    gt = genotype_from_gt(gt)
    if gt in ("", ".", "./.", ".|."):
        return "Missing"
    mapping = {"0": ref, "1": alt}
    parts = gt.split("/")
    if len(parts) != 2:
        return gt
    return f"{mapping.get(parts[0], parts[0])}/{mapping.get(parts[1], parts[1])}"


def calculate_file_hash(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_vcf(file_bytes: bytes) -> Tuple[pd.DataFrame, List[str]]:
    warnings = []
    rows = []
    text = file_bytes.decode("utf-8", errors="replace")

    header_samples = []
    for line in text.splitlines():
        if not line or line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            cols = line.lstrip("#").split("\t")
            header_samples = cols[9:]
            continue

        if line.startswith("#"):
            continue

        parts = line.split("\t")
        if len(parts) < 8:
            warnings.append("Skipped malformed VCF row with fewer than 8 columns.")
            continue

        chrom, pos, rsid, ref, alt, qual, flt, info = parts[:8]
        fmt = parts[8] if len(parts) > 8 else ""
        sample_values = parts[9:] if len(parts) > 9 else []

        sample_name = header_samples[0] if header_samples else "Sample"
        sample_value = sample_values[0] if sample_values else ""

        format_keys = fmt.split(":") if fmt else []
        format_map = {}
        if sample_value and format_keys:
            vals = sample_value.split(":")
            format_map = dict(zip(format_keys, vals))

        gt = format_map.get("GT", "")
        dp = format_map.get("DP", "")
        gq = format_map.get("GQ", "")

        for alt_allele in alt.split(","):
            rid = normalize_rs(rsid)
            ann = VARIANT_MAP.get(rid, {})
            gene = ann.get("gene", "")
            if not gene:
                # Try gene symbols in INFO if present.
                upper_info = info.upper()
                for g in GENES:
                    if g in upper_info:
                        gene = g
                        break

            rows.append({
                "CHROM": chrom,
                "POS": pos,
                "ID": rsid,
                "REF": ref,
                "ALT": alt_allele,
                "QUAL": qual,
                "FILTER": flt,
                "GT": gt,
                "Genotype": allele_pair_from_gt(gt, ref, alt_allele),
                "DP": dp,
                "GQ": gq,
                "Gene": gene,
                "Known interpretation": ann.get("label", ""),
                "Allele mapping": ann.get("allele", ""),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        warnings.append("No variant records were parsed from the VCF.")
    if not header_samples:
        warnings.append("No sample column was detected; variant-level parsing only.")
    return df, warnings


def parse_csv_txt(uploaded_file) -> Tuple[pd.DataFrame, List[str]]:
    warnings = []
    raw = uploaded_file.getvalue()
    name = uploaded_file.name.lower()

    try:
        if name.endswith(".csv"):
            df = pd.read_csv(io.BytesIO(raw))
        else:
            text = raw.decode("utf-8", errors="replace")
            # Try comma, tab, then whitespace.
            try:
                df = pd.read_csv(io.StringIO(text), sep=None, engine="python")
            except Exception:
                df = pd.read_csv(io.StringIO(text), sep="\t")
    except Exception as exc:
        return pd.DataFrame(), [f"Could not parse file: {exc}"]

    df.columns = [str(c).strip() for c in df.columns]

    # Normalize common column names.
    aliases = {}
    for c in df.columns:
        lc = c.lower().replace("_", "").replace("-", "")
        if lc in ("gene", "genesymbol", "symbol"):
            aliases[c] = "Gene"
        elif lc in ("rsid", "rs", "dbsnp"):
            aliases[c] = "rsID"
        elif lc in ("genotype", "gt", "call"):
            aliases[c] = "GT"
        elif lc in ("chrom", "chr", "chromosome"):
            aliases[c] = "CHROM"
        elif lc in ("pos", "position"):
            aliases[c] = "POS"
        elif lc == "ref":
            aliases[c] = "REF"
        elif lc == "alt":
            aliases[c] = "ALT"
        elif lc == "qual":
            aliases[c] = "QUAL"
        elif lc == "filter":
            aliases[c] = "FILTER"
    df = df.rename(columns=aliases)

    if "Gene" in df.columns:
        df["Gene"] = df["Gene"].map(normalize_gene)
    else:
        df["Gene"] = ""

    if "rsID" in df.columns:
        df["rsID"] = df["rsID"].map(normalize_rs)
    else:
        df["rsID"] = ""

    if "GT" in df.columns:
        df["GT"] = df["GT"].map(genotype_from_gt)
    else:
        df["GT"] = ""

    labels, allele_maps = [], []
    for _, row in df.iterrows():
        ann = VARIANT_MAP.get(row.get("rsID", ""), {})
        labels.append(ann.get("label", ""))
        allele_maps.append(ann.get("allele", ""))
    df["Known interpretation"] = labels
    df["Allele mapping"] = allele_maps

    if df.empty:
        warnings.append("The uploaded table contains no rows.")
    return df, warnings


def qc_variant_table(df: pd.DataFrame) -> Dict[str, object]:
    qc = {
        "status": "PASS",
        "variant_count": int(len(df)),
        "recognized_variants": 0,
        "missing_genotypes": 0,
        "low_quality": 0,
        "warnings": [],
    }

    if df.empty:
        qc["status"] = "FAIL"
        qc["warnings"].append("No variants available for analysis.")
        return qc

    if "Known interpretation" in df.columns:
        qc["recognized_variants"] = int((df["Known interpretation"].fillna("") != "").sum())

    if "GT" in df.columns:
        qc["missing_genotypes"] = int(
            df["GT"].fillna("").astype(str).isin(["", ".", "./.", ".|."]).sum()
        )

    if "QUAL" in df.columns:
        numeric = pd.to_numeric(df["QUAL"], errors="coerce")
        qc["low_quality"] = int((numeric.notna() & (numeric < 20)).sum())

    if qc["missing_genotypes"] > 0:
        qc["warnings"].append(
            f"{qc['missing_genotypes']} variant record(s) have missing genotype calls."
        )
    if qc["low_quality"] > 0:
        qc["warnings"].append(
            f"{qc['low_quality']} variant record(s) have QUAL < 20; review before interpretation."
        )
    if qc["recognized_variants"] == 0:
        qc["warnings"].append(
            "No variants matched the built-in evidence map. This does not mean the sample is negative."
        )

    if qc["warnings"]:
        qc["status"] = "WARNING"
    return qc


# -----------------------------
# Phenotype engine
# -----------------------------
def phenotype_from_diplotype(gene: str, diplotype: str) -> Tuple[str, str]:
    d = diplotype.replace(" ", "").upper()
    gene = normalize_gene(gene)

    if gene == "CYP2D6":
        if "*1XN" in d or "XN" in d:
            return "Ultrarapid metabolizer", "Manual/validated copy-number-aware diplotype input"
        if d in {"*4/*4", "*4/*5", "*5/*4", "*5/*5"}:
            return "Poor metabolizer", "Two no/low-function alleles reported"
        if d in {"*1/*1", "*1/*2", "*2/*2"}:
            return "Normal metabolizer", "Normal-function diplotype supplied"
        if "*4" in d or "*10" in d:
            return "Intermediate metabolizer", "Reduced/no-function allele combination requires validated allele activity assignment"
        return "Unresolved", "Diplotype not recognized by the conservative built-in rule set"

    if gene == "CYP2C19":
        if d == "*2/*2":
            return "Poor metabolizer", "Two loss-of-function alleles reported"
        if d in {"*1/*2", "*2/*17"}:
            return "Intermediate metabolizer", "Reduced-function allele present"
        if d in {"*1/*1"}:
            return "Normal metabolizer", "Normal-function diplotype supplied"
        if "*17" in d:
            return "Rapid metabolizer", "Increased-function allele reported; confirm exact diplotype and guideline applicability"
        return "Unresolved", "Diplotype not recognized by the conservative built-in rule set"

    if gene == "CYP3A4":
        if "*22" in d:
            return "Reduced function", "CYP3A4*22-associated finding; drug-specific interpretation requ
