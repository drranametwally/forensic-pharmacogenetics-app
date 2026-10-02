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
# R-FINGERPRINT DNA SUITE[span_0](start_span)[span_0](end_span)
# FORENSIC PHARMACOGENETICS & TOXICOLOGY DECISION SUPPORT
# Version 3.3 — Complete with Caching, Pre-validation & Enhanced Toxicity/Addiction Markers
# ============================================================

st.set_page_config(
    page_title="R-Fingerprint DNA Suite",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

APP_VERSION = "3.3"
BUILD_DATE = datetime.now().strftime("%Y-%m-%d")

GENES = [
    "CYP2D6", "CYP2C19", "CYP3A4", "OPRM1", "ABCB1",
    "CYP2B6", "SLCO1B1", "HLA-B", "DRD2", "HTR2A", "CYP2C9", "COMT"
]

# قاموس المتغيرات الموسع ليشمل استقلاب الأدوية والسموم والإدمان
VARIANT_MAP = {
    # الأدوية الأساسية واستقلاب السموم العصبية
    "rs3892097": {"gene": "CYP2D6", "label": "CYP2D6*4-defining variant", "allele": "*4"},
    "rs35742686": {"gene": "CYP2D6", "label": "CYP2D6*3-defining variant", "allele": "*3"},
    "rs1065852": {"gene": "CYP2D6", "label": "CYP2D6*10-associated variant", "allele": "*10"},
    "rs5030865": {"gene": "CYP2D6", "label": "CYP2D6*14-associated variant", "allele": "*14"},
    "rs1799971": {"gene": "OPRM1", "label": "OPRM1 A118G (Opioid response & receptor sensitivity)", "allele": "A118G"},
    "rs4244285": {"gene": "CYP2C19", "label": "CYP2C19*2", "allele": "*2"},
    "rs12248560": {"gene": "CYP2C19", "label": "CYP2C19*17-associated variant", "allele": "*17"},
    "rs28399504": {"gene": "CYP2C19", "label": "CYP2C19*2-associated variant", "allele": "*2"},
    "rs35599367": {"gene": "CYP3A4", "label": "CYP3A4*22-associated variant", "allele": "*22"},
    "rs4149056": {"gene": "SLCO1B1", "label": "SLCO1B1 c.521T>C", "allele": "*5/*15-associated"},
    "rs1045642": {"gene": "ABCB1", "label": "ABCB1 3435C>T", "allele": "3435C>T"},
    "rs2395029": {"gene": "HLA-B", "label": "HLA-linked marker; requires validated HLA interpretation", "allele": "marker"},
    
    # الـ 5 متغيرات الجديدة لربط البيانات بالاستبيان الميداني وحالات التسمم والإدمان
    "rs1800497": {"gene": "DRD2", "label": "DRD2/ANKK1 Taq1A (Addiction/Dopamine response)", "allele": "A1/A2"},
    "rs6311": {"gene": "HTR2A", "label": "HTR2A -1038 C>T (Serotonergic/Toxicity response)", "allele": "T"},
    "rs2032582": {"gene": "ABCB1", "label": "ABCB1 c.2677G>T/A (Blood-brain barrier transport)", "allele": "T/A"},
    "rs1799853": {"gene": "CYP2C9", "label": "CYP2C9*2 (Metabolism of co-ingested drugs)", "allele": "*2"},
    "rs4680": {"gene": "COMT", "label": "COMT Val158Met (Pain threshold & neurotransmitter clearance)", "allele": "Met"},
}

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
            "mechanism": "Reduced conversion of codeine to morphine may reduce analgesic effect.",
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
            "mechanism": "CYP2D6 phenotype may alter systemic concentrations.",
            "evidence": "FDA pharmacogenetic association table",
            "source": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
        }
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
            "mechanism": "CYP2B6 variation can influence methadone disposition.",
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
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-title">🧬 R-Fingerprint DNA Suite</div>', unsafe_allow_html=True)
st.markdown('<div class="subtitle">Forensic Pharmacogenetics & Toxicology Decision-Support — Master Thesis Edition (v3.3)</div>', unsafe_allow_html=True)

# -----------------------------
# Utility functions & Cached Parsers
# -----------------------------
def normalize_gene(value: str) -> str:
    value = str(value).strip().upper().replace(" ", "")
    aliases = {
        "CYP2D6": "CYP2D6", "CYP2C19": "CYP2C19", "CYP3A4": "CYP3A4",
        "OPRM1": "OPRM1", "ABCB1": "ABCB1", "CYP2B6": "CYP2B6",
        "SLCO1B1": "SLCO1B1", "HLA-B": "HLA-B", "HLAB": "HLA-B",
        "DRD2": "DRD2", "HTR2A": "HTR2A", "CYP2C9": "CYP2C9", "COMT": "COMT"
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


# 1. Caching & Pre-validation for VCF
@st.cache_data
def parse_vcf_cached(file_bytes: bytes) -> Tuple[pd.DataFrame, List[str]]:
    warnings = []
    rows = []
    text = file_bytes.decode("utf-8", errors="replace")

    if "#CHROM" not in text and not any(line.startswith("##fileformat") for line in text.splitlines()[:5]):
        return pd.DataFrame(), ["Error: The uploaded file does not appear to be a valid VCF format."]

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
        sample_value = parts[9] if len(parts) > 9 else ""

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
                upper_info = info.upper()
                for g in GENES:
                    if g in upper_info:
                        gene = g
                        break

            rows.append({
                "CHROM": chrom, "POS": pos, "ID": rsid, "REF": ref,
                "ALT": alt_allele, "QUAL": qual, "FILTER": flt, "GT": gt,
                "Genotype": allele_pair_from_gt(gt, ref, alt_allele),
                "DP": dp, "GQ": gq, "Gene": gene,
                "Known interpretation": ann.get("label", ""),
                "Allele mapping": ann.get("allele", ""),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        warnings.append("No variant records were parsed from the VCF.")
    return df, warnings


# 1. Caching & Pre-validation for CSV/TXT
@st.cache_data
def parse_csv_txt_cached(raw: bytes, filename: str) -> Tuple[pd.DataFrame, List[str]]:
    warnings = []
    try:
        if filename.endswith(".csv"):
            df = pd.read_csv(io.BytesIO(raw))
        else:
            text = raw.decode("utf-8", errors="replace")
            try:
                df = pd.read_csv(io.StringIO(text), sep=None, engine="python")
            except Exception:
                df = pd.read_csv(io.StringIO(text), sep="\t")
    except Exception as exc:
        return pd.DataFrame(), [f"Could not parse file: {exc}"]

    if df.empty:
        return pd.DataFrame(), ["Error: The uploaded table contains no rows."]

    df.columns = [str(c).strip() for c in df.columns]

    aliases = {}
    for c in df.columns:
        lc = c.lower().replace("_", "").replace("-", "")
        if lc in ("gene", "genesymbol", "symbol"): aliases[c] = "Gene"
        elif lc in ("rsid", "rs", "dbsnp"): aliases[c] = "rsID"
        elif lc in ("genotype", "gt", "call"): aliases[c] = "GT"
        elif lc in ("chrom", "chr", "chromosome"): aliases[c] = "CHROM"
        elif lc in ("pos", "position"): aliases[c] = "POS"

    df = df.rename(columns=aliases)

    if "rsID" not in df.columns and "Gene" not in df.columns:
        warnings.append("Warning: Missing standard columns ('rsID' or 'Gene'). Interpretation might be limited.")

    df["Gene"] = df["Gene"].map(normalize_gene) if "Gene" in df.columns else ""
    df["rsID"] = df["rsID"].map(normalize_rs) if "rsID" in df.columns else ""
    df["GT"] = df["GT"].map(genotype_from_gt) if "GT" in df.columns else ""

    labels, allele_maps = [], []
    for _, row in df.iterrows():
        ann = VARIANT_MAP.get(row.get("rsID", ""), {})
        labels.append(ann.get("label", ""))
        allele_maps.append(ann.get("allele", ""))
    df["Known interpretation"] = labels
    df["Allele mapping"] = allele_maps

    return df, warnings


def qc_variant_table(df: pd.DataFrame) -> Dict[str, object]:
    qc = {
        "status": "PASS", "variant_count": int(len(df)),
        "recognized_variants": 0, "missing_genotypes": 0,
        "low_quality": 0, "warnings": [],
    }
    if df.empty:
        qc["status"] = "FAIL"
        qc["warnings"].append("No variants available for analysis.")
        return qc

    if "Known interpretation" in df.columns:
        qc["recognized_variants"] = int((df["Known interpretation"].fillna("") != "").sum())
    if "GT" in df.columns:
        qc["missing_genotypes"] = int(df["GT"].fillna("").astype(str).isin(["", ".", "./.", ".|."]).sum())
    if "QUAL" in df.columns:
        numeric = pd.to_numeric(df["QUAL"], errors="coerce")
        qc["low_quality"] = int((numeric.notna() & (numeric < 20)).sum())

    if qc["missing_genotypes"] > 0:
        qc["warnings"].append(f"{qc['missing_genotypes']} variant record(s) have missing genotype calls.")
    if qc["low_quality"] > 0:
        qc["warnings"].append(f"{qc['low_quality']} variant record(s) have low quality score (QUAL < 20).")
    
    if qc["warnings"]:
        qc["status"] = "WARNING"
    return qc


# -----------------------------
# Streamlit UI Layout
# -----------------------------
with st.sidebar:
    st.header("⚙️ Case Setup")
    case_id = st.text_input("Case ID", value=f"CASE-{datetime.now().strftime('%Y%m%d-%H%M%S')}")
    drug = st.selectbox("Drug assessed", ["Codeine", "Tramadol", "Clopidogrel", "Amitriptyline", "Omeprazole", "Methadone", "Other / not in evidence map"])
    input_method = st.radio("Input method", ["VCF / genomic file", "CSV / TXT table", "Interactive panel"])

st.markdown('<div class="section-title">1. Genomic Input & QC (Cached with Toxicity & Addiction Panel)</div>', unsafe_allow_html=True)

uploaded = None
variant_df = pd.DataFrame()
parse_warnings = []
file_hash = ""

if input_method == "VCF / genomic file":
    uploaded = st.file_uploader("Upload VCF", type=["vcf", "txt"])
elif input_method == "CSV / TXT table":
    uploaded = st.file_uploader("Upload CSV / TXT", type=["csv", "txt"])

if uploaded is not None:
    raw = uploaded.getvalue()
    file_hash = calculate_file_hash(raw)

    if input_method == "VCF / genomic file":
        variant_df, parse_warnings = parse_vcf_cached(raw)
    else:
        variant_df, parse_warnings = parse_csv_txt_cached(raw, uploaded.name.lower())

    qc = qc_variant_table(variant_df)
else:
    qc = {"status": "NOT RUN", "variant_count": 0, "recognized_variants": 0, "missing_genotypes": 0, "warnings": ["No genomic file supplied."]}

qc_cols = st.columns(4)
qc_cols[0].metric("QC", qc["status"])
qc_cols[1].metric("Variants parsed", qc["variant_count"])
qc_cols[2].metric("Evidence matches", qc["recognized_variants"])
qc_cols[3].metric("Missing GT", qc["missing_genotypes"])

if parse_warnings:
    for w in parse_warnings:
        if "Error" in w:
            st.error(w)
        else:
            st.warning(w)

if not variant_df.empty:
    with st.expander("Parsed genomic records (Including Addiction & Forensic Markers)", expanded=False):
        st.dataframe(variant_df, use_container_width=True, hide_index=True)

