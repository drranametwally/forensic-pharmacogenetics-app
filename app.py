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


def normalize_diplotype(value: str) -> str:
    """Normalize common star-allele/diplotype spellings without inventing a call."""
    if value is None:
        return ""
    d = str(value).strip().upper()
    if d in {"", "NAN", "NONE", "N/A", "NA", "."}:
        return ""
    d = d.replace("\\", "/").replace("|", "/")
    d = re.sub(r"\s+", "", d)
    return d


def parse_star_diplotype(diplotype: str) -> Optional[List[str]]:
    """Return two star-allele tokens only when the supplied value is structurally parseable."""
    d = normalize_diplotype(diplotype)
    if not d:
        return None
    parts = d.split("/")
    if len(parts) != 2 or not all(re.fullmatch(r"\*[^/]+", p) for p in parts):
        return None
    return parts


def allele_count_in_diplotype(diplotype: str, allele: str) -> Optional[int]:
    parts = parse_star_diplotype(diplotype)
    if parts is None:
        return None
    return sum(p == allele.upper() for p in parts)


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
        elif lc in ("diplotype", "diplotypes", "stardiplotype", "staralleles", "starallele"):
            aliases[c] = "Diplotype"
    df = df.rename(columns=aliases)

    # Preserve a supplied diplotype exactly as a normalized input field.
    # It is evidence supplied by the laboratory/caller, not something inferred here.
    if "Diplotype" in df.columns:
        df["Diplotype"] = df["Diplotype"].map(normalize_diplotype)

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


def validate_genotype_diplotype_consistency(df: pd.DataFrame) -> List[Dict[str, str]]:
    """Check only defensible genotype↔diplotype relationships; never infer missing alleles."""
    conflicts = []
    if df.empty or "Diplotype" not in df.columns:
        return conflicts

    # These markers are treated as strict consistency checks because the mapped
    # variant corresponds directly to the named star allele in the supplied diplotype.
    strict_map = {
        "rs3892097": "*4",
        "rs35742686": "*3",
        "rs4244285": "*2",
        "rs1799853": "*2",
        "rs1057910": "*3",
        "rs35599367": "*22",
    }

    for gene, sub in df.groupby(df["Gene"].astype(str).str.upper()):
        values = [v for v in sub["Diplotype"].tolist() if normalize_diplotype(v)]
        if not values:
            continue
        unique_values = list(dict.fromkeys(normalize_diplotype(v) for v in values))
        if len(unique_values) > 1:
            conflicts.append({
                "gene": gene,
                "type": "Multiple diplotypes",
                "message": f"Conflicting diplotypes supplied for {gene}: {', '.join(unique_values)}.",
            })
            continue

        diplotype = unique_values[0]
        if parse_star_diplotype(diplotype) is None and gene in {"CYP2D6", "CYP2C19", "CYP3A4", "CYP2B6", "CYP2C9"}:
            conflicts.append({
                "gene": gene,
                "type": "Unparseable diplotype",
                "message": f"Supplied diplotype '{diplotype}' could not be structurally validated; phenotype assignment must be reviewed.",
            })
            continue

        for _, row in sub.iterrows():
            rsid = str(row.get("rsID", "")).strip().lower()
            expected_allele = strict_map.get(rsid)
            if not expected_allele:
                continue
            gt = genotype_from_gt(row.get("GT", ""))
            if gt not in {"0/0", "0/1", "1/1"}:
                continue
            copies = allele_count_in_diplotype(diplotype, expected_allele)
            if copies is None:
                continue
            expected_gt = {0: "0/0", 1: "0/1", 2: "1/1"}[copies]
            if gt != expected_gt:
                conflicts.append({
                    "gene": gene,
                    "type": "GT/Diplotype conflict",
                    "message": (
                        f"{gene} {rsid}: genotype {gt} conflicts with reported diplotype {diplotype}; "
                        f"the consistency check expects {expected_gt} for {expected_allele}."
                    ),
                })
    return conflicts


def qc_variant_table(df: pd.DataFrame) -> Dict[str, object]:
    qc = {
        "status": "PASS",
        "interpretation_status": "REVIEW REQUIRED",
        "variant_count": int(len(df)),
        "recognized_variants": 0,
        "positive_mapped_variants": 0,
        "missing_genotypes": 0,
        "low_quality": 0,
        "invalid_genotypes": 0,
        "diplotype_conflicts": [],
        "warnings": [],
    }

    if df.empty:
        qc["status"] = "FAIL"
        qc["interpretation_status"] = "NOT INTERPRETABLE"
        qc["warnings"].append("No variants available for analysis.")
        return qc

    if "Known interpretation" in df.columns:
        mapped = df["Known interpretation"].fillna("") != ""
        qc["recognized_variants"] = int(mapped.sum())
        if "GT" in df.columns:
            gt_for_map = df["GT"].fillna("").astype(str).map(genotype_from_gt)
            qc["positive_mapped_variants"] = int((mapped & gt_for_map.isin(["0/1", "1/1"])).sum())

    if "GT" in df.columns:
        gt_series = df["GT"].fillna("").astype(str).map(genotype_from_gt)
        qc["missing_genotypes"] = int(gt_series.isin(["", ".", "./.", ".|."]).sum())
        qc["invalid_genotypes"] = int((~gt_series.isin(["", ".", "./.", ".|.", "0/0", "0/1", "1/1"])).sum())

    if "QUAL" in df.columns:
        numeric = pd.to_numeric(df["QUAL"], errors="coerce")
        qc["low_quality"] = int((numeric.notna() & (numeric < 20)).sum())

    qc["diplotype_conflicts"] = validate_genotype_diplotype_consistency(df)

    if qc["missing_genotypes"] > 0:
        qc["warnings"].append(
            f"{qc['missing_genotypes']} variant record(s) have missing genotype calls."
        )
    if qc["invalid_genotypes"] > 0:
        qc["warnings"].append(
            f"{qc['invalid_genotypes']} genotype call(s) are unsupported/non-biallelic; review before interpretation."
        )
    if qc["low_quality"] > 0:
        qc["warnings"].append(
            f"{qc['low_quality']} variant record(s) have QUAL < 20; review before interpretation."
        )
    if qc["recognized_variants"] == 0:
        qc["warnings"].append(
            "No variants matched the built-in evidence map. This does not mean the sample is negative."
        )
    if qc["diplotype_conflicts"]:
        qc["warnings"].append(
            f"{len(qc['diplotype_conflicts'])} genotype/diplotype consistency issue(s) detected. Automated phenotype interpretation is blocked for affected genes."
        )

    if qc["warnings"]:
        qc["status"] = "WARNING"

    # PASS means the file passed basic structural QC only. It is deliberately
    # not a statement of analytical validity or clinical validity.
    if qc["missing_genotypes"] == 0 and qc["invalid_genotypes"] == 0 and not qc["diplotype_conflicts"]:
        qc["interpretation_status"] = "STRUCTURALLY CONSISTENT — PROFESSIONAL REVIEW REQUIRED"
    else:
        qc["interpretation_status"] = "REVIEW REQUIRED — INTERPRETATION MAY BE BLOCKED"
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
            return "Reduced function", "CYP3A4*22-associated finding; drug-specific interpretation required"
        if "*1" in d:
            return "No actionable reduced-function allele detected", "Limited to supplied diplotype"
        return "Unresolved", "Drug-specific interpretation required"

    if gene == "CYP2B6":
        if "*6" in d or "*18" in d:
            return "Reduced function", "Common reduced-function allele reported; exact drug-specific interpretation required"
        if "*1/*1" in d:
            return "Normal function", "Normal-function diplotype supplied"
        return "Unresolved", "Diplotype not recognized"

    if gene == "CYP2C9":
        if d in {"*2/*2", "*2/*3", "*3/*2", "*3/*3"}:
            return "Poor / markedly reduced function", "Two reduced-function alleles reported; drug-specific interpretation required"
        if d in {"*1/*2", "*2/*1", "*1/*3", "*3/*1"}:
            return "Intermediate / reduced function", "One reduced-function allele reported; drug-specific interpretation required"
        if d == "*1/*1":
            return "Normal metabolizer", "Normal-function diplotype supplied"
        return "Unresolved", "CYP2C9 diplotype not recognized"

    if gene in {"OPRM1", "ABCB1", "SLCO1B1", "HLA-B", "HTR2A", "ANKK1", "COMT"}:
        return "Variant-level finding", "A general phenotype should not be inferred without a validated gene-specific interpretation model"

    return "Unresolved", "No built-in phenotype model"


def variant_level_findings(df: pd.DataFrame) -> List[Dict[str, str]]:
    findings = []
    if df.empty:
        return findings

    for _, row in df.iterrows():
        rsid = str(row.get("rsID", "")).strip()
        ann = VARIANT_MAP.get(rsid)
        if not ann:
            continue

        gt = genotype_from_gt(row.get("GT", ""))
        if gt == "0/0":
            finding = "Reference genotype; mapped alternate allele not detected"
        elif gt in {"0/1", "1/1"}:
            finding = ann["label"]
        else:
            finding = "Genotype present but not interpretable by the current biallelic rule set"
        findings.append({
            "gene": ann["gene"],
            "rsID": rsid,
            "genotype": gt,
            "finding": finding,
            "allele_mapping": ann["allele"],
        })
    return findings


def receptor_assessment(findings: List[Dict[str, str]], manual_oprm1: str) -> Tuple[str, str]:
    if manual_oprm1 == "A118G variant reported":
        return (
            "Variant-level receptor finding",
            "OPRM1 A118G is present in the supplied interpretation. Drug-specific clinical effect must be confirmed; presence alone is not a validated universal overdose predictor."
        )

    for f in findings:
        if f["gene"] == "OPRM1" and f["rsID"] == "rs1799971":
            gt = str(f.get("genotype", "")).strip()

            if gt in {"0/1", "1/1"}:
                return (
                    "Variant-level receptor finding",
                    "OPRM1 A118G detected. Do not convert this finding into a universal 'high sensitivity' or overdose label without drug-specific validated evidence."
                )

            if gt == "0/0":
                return (
                    "No actionable receptor finding",
                    "OPRM1 rs1799971 is reference (0/0); the mapped alternate allele was not detected."
                )

            return (
                "Receptor finding requires review",
                "OPRM1 rs1799971 was detected, but the genotype is not interpretable by the current biallelic rule set."
            )

    return (
        "No actionable receptor finding",
        "No supported OPRM1 finding was identified in the current input."
    )


def metabolic_summary(phenotypes: Dict[str, Tuple[str, str]]) -> Tuple[str, str]:
    if not phenotypes:
        return "Insufficient genomic interpretation", "No gene-level phenotype was established."

    values = [p[0] for p in phenotypes.values()]
    if "Poor metabolizer" in values:
        return "Compromised / Poor metabolism identified", "At least one gene has a poor-metabolizer phenotype."
    if "Intermediate metabolizer" in values or "Reduced function" in values:
        return "Reduced / Intermediate metabolic function", "At least one gene has reduced metabolic function."
    if "Ultrarapid metabolizer" in values or "Rapid metabolizer" in values:
        return "Increased metabolic function identified", "At least one gene has an increased-function phenotype."
    if "Normal metabolizer" in values or "Normal function" in values:
        return "No reduced function identified in interpreted genes", "Normal function was identified for the supplied diplotype(s)."
    return "Variant findings require review", "No validated gene-level phenotype was established."


def risk_engine(drug: str, phenotypes: Dict[str, Tuple[str, str]], receptor_text: str,
                exposure_factors: Dict[str, object]) -> Tuple[str, List[Dict[str, str]], str]:
    matches = []

    for rule in DRUG_RULES.get(drug, []):
        gene = rule["gene"]
        pheno = phenotypes.get(gene, ("", ""))[0]
        if pheno in rule["phenotypes"]:
            matches.append({
                "gene": gene,
                "indicator": rule["indicator"],
                "severity": rule["severity"],
                "mechanism": rule["mechanism"],
                "evidence": rule["evidence"],
                "source": rule["source"],
            })

    # Do NOT convert exposure history alone into genetic risk.
    # It is used as contextual information.
    if not matches:
        return (
            "No drug-specific actionable PGx signal established",
            [],
            "No matching rule was found for the selected drug and interpreted phenotype. This is not evidence of zero toxicological risk."
        )

    max_rank = max(SEVERITY_RANK[m["severity"]] for m in matches)
    if max_rank >= 3:
        overall = "High pharmacogenetic concern — professional review required"
    elif max_rank == 2:
        overall = "Moderate pharmacogenetic concern — professional review required"
    else:
        overall = "Low pharmacogenetic concern"

    context = (
        f"Genetic interpretation is combined with clinical/toxicological context. "
        f"Reported overdose history: {exposure_factors['overdose_history']}. "
        f"Co-medication/polypharmacy: {exposure_factors['polypharmacy']}."
    )
    return overall, matches, context


def manual_panel_interpretation(values: Dict[str, str]) -> Tuple[Dict[str, Tuple[str, str]], List[Dict[str, str]]]:
    phenotypes = {}
    findings = []

    for gene, value in values.items():
        if value == "Not assessed":
            continue

        if gene in {"CYP2D6", "CYP2C19", "CYP3A4", "CYP2B6", "CYP2C9"}:
            pheno, reason = phenotype_from_diplotype(gene, value)
            phenotypes[gene] = (pheno, reason)
            findings.append({
                "gene": gene,
                "rsID": "Manual",
                "genotype": value,
                "finding": pheno,
                "allele_mapping": reason,
            })
        elif gene == "OPRM1" and value == "A118G variant reported":
            findings.append({
                "gene": "OPRM1",
                "rsID": "Manual",
                "genotype": value,
                "finding": "A118G variant-level finding",
                "allele_mapping": "Drug-specific interpretation required",
            })
        elif gene == "ABCB1" and value != "Not assessed":
            for rsid in (["rs2032582", "rs1045642"] if value == "Both variants reported" else (["rs2032582"] if "rs2032582" in value else ["rs1045642"])):
                ann = VARIANT_MAP[rsid]
                findings.append({"gene": "ABCB1", "rsID": rsid, "genotype": value, "finding": ann["label"], "allele_mapping": ann["allele"]})
        elif gene in {"HTR2A", "ANKK1", "COMT"} and value != "Not assessed":
            mapping = {
                "HTR2A": {"rs6311 reported": ["rs6311"], "rs6313 reported": ["rs6313"], "Both variants reported": ["rs6311", "rs6313"]},
                "ANKK1": {"rs1800497 reported": ["rs1800497"]},
                "COMT": {"rs4680 Val158Met reported": ["rs4680"]},
            }
            for rsid in mapping.get(gene, {}).get(value, []):
                ann = VARIANT_MAP[rsid]
                findings.append({"gene": gene, "rsID": rsid, "genotype": value, "finding": ann["label"], "allele_mapping": ann["allele"]})
    return phenotypes, findings


# -----------------------------
# PDF report
# -----------------------------
def make_pdf(report: Dict) -> bytes:
    if not REPORTLAB_OK:
        raise RuntimeError("reportlab is not installed. Add reportlab to requirements.txt.")

    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        rightMargin=14 * mm,
        leftMargin=14 * mm,
        topMargin=14 * mm,
        bottomMargin=14 * mm,
    )

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="CenterTitle", parent=styles["Title"], alignment=TA_CENTER, fontSize=17, leading=21
    ))
    styles.add(ParagraphStyle(
        name="Small", parent=styles["BodyText"], fontSize=8, leading=10
    ))

    story = []
    story.append(Paragraph("R-Fingerprint DNA Suite", styles["CenterTitle"]))
    story.append(Paragraph("Forensic Pharmacogenetics & Toxicology Decision-Support Report", styles["Heading2"]))
    story.append(Spacer(1, 6))

    meta = [
        ["Case ID", report["case_id"]],
        ["Case type", report["case_type"]],
        ["Selected drug/substance", report["drug"]],
        ["Generated", report["generated_at"]],
        ["Software version", APP_VERSION],
        ["Input SHA-256", report["file_hash"] or "No genomic file supplied"],
    ]
    t = Table(meta, colWidths=[45 * mm, 135 * mm])
    t.setStyle(TableStyle([
        ("GRID", (0,0), (-1,-1), 0.4, colors.grey),
        ("BACKGROUND", (0,0), (0,-1), colors.whitesmoke),
        ("VALIGN", (0,0), (-1,-1), "TOP"),
    ]))
    story.append(t)
    story.append(Spacer(1, 10))

    story.append(Paragraph("1. QC Status", styles["Heading2"]))
    story.append(Paragraph(
        f"Structural QC: <b>{report['qc']['status']}</b> | Variants parsed: {report['qc']['variant_count']} | "
        f"Mapped loci: {report['qc']['recognized_variants']} | Positive mapped: {report['qc'].get('positive_mapped_variants', 0)} | "
        f"Missing GT: {report['qc']['missing_genotypes']} | Invalid GT: {report['qc'].get('invalid_genotypes', 0)}",
        styles["BodyText"]
    ))
    story.append(Paragraph(
        f"Interpretation gate: <b>{report['qc'].get('interpretation_status', 'REVIEW REQUIRED')}</b>. "
        "Structural QC status is not analytical validation or clinical validation.",
        styles["Small"]
    ))
    for warning in report["qc"]["warnings"]:
        story.append(Paragraph(f"• {warning}", styles["Small"]))
    for conflict in report["qc"].get("diplotype_conflicts", []):
        story.append(Paragraph(f"• {conflict['message']}", styles["Small"]))

    story.append(Spacer(1, 8))
    story.append(Paragraph("2. Primary Results", styles["Heading2"]))
    primary = [
        ["Metabolic status", report["metabolic_status"]],
        ["Receptor assessment", report["receptor_status"]],
        ["Toxicology / PGx indicator", report["risk_status"]],
    ]
    t = Table(primary, colWidths=[55 * mm, 125 * mm])
    t.setStyle(TableStyle([
        ("GRID", (0,0), (-1,-1), 0.4, colors.grey),
        ("BACKGROUND", (0,0), (0,-1), colors.whitesmoke),
        ("VALIGN", (0,0), (-1,-1), "TOP"),
    ]))
    story.append(t)
    story.append(Spacer(1, 8))

    story.append(Paragraph("3. Gene Findings & Interpretation Basis", styles["Heading2"]))
    rows = [["Gene", "Variant / Genotype", "Reported diplotype", "Phenotype / interpretation"]]
    for f in report["findings"]:
        gene = f["gene"]
        pheno = report.get("phenotypes", {}).get(gene, ("Not assigned", ""))
        reported_diplotype = "—"
        if gene in {"CYP2D6", "CYP2C19", "CYP3A4", "CYP2B6", "CYP2C9"}:
            # Find the explicit reported diplotype in the finding's source text when available.
            reason = str(pheno[1]) if isinstance(pheno, tuple) and len(pheno) > 1 else ""
            m = re.search(r"Reported diplotype ([^;]+)", reason)
            if m:
                reported_diplotype = m.group(1)
        interpretation = f["finding"]
        if pheno[0] != "Not assigned":
            interpretation += f"; Phenotype: {pheno[0]}"
        rows.append([
            gene,
            f"{f['rsID']} | {f['genotype']}",
            reported_diplotype,
            interpretation,
        ])
    if len(rows) == 1:
        rows.append(["—", "—", "No supported finding available"])
    t = Table(rows, colWidths=[25 * mm, 48 * mm, 38 * mm, 64 * mm])
    t.setStyle(TableStyle([
        ("GRID", (0,0), (-1,-1), 0.4, colors.grey),
        ("BACKGROUND", (0,0), (-1,0), colors.whitesmoke),
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("FONTSIZE", (0,0), (-1,-1), 8),
    ]))
    story.append(t)
    story.append(Spacer(1, 8))

    story.append(Paragraph("4. Drug-Specific Evidence", styles["Heading2"]))
    for m in report["matches"]:
        story.append(Paragraph(
            f"<b>{m['gene']} — {m['indicator']} ({m['severity']})</b><br/>"
            f"{m['mechanism']}<br/>"
            f"Evidence: {m['evidence']}<br/>"
            f"Source: {m['source']}",
            styles["Small"]
        ))
        story.append(Spacer(1, 4))
    if not report["matches"]:
        story.append(Paragraph(
            "No drug-specific actionable pharmacogenetic rule matched the interpreted phenotype.",
            styles["BodyText"]
        ))

    story.append(Spacer(1, 8))
    story.append(Paragraph("5. Clinical/Toxicological Context", styles["Heading2"]))
    story.append(Paragraph(report["context"], styles["BodyText"]))

    story.append(Spacer(1, 8))
    story.append(Paragraph("6. Limitations and Required Review", styles["Heading2"]))
    limitations = [
        "This report is decision-support output, not an autonomous diagnosis or prescribing instruction.",
        "A single variant is not automatically equivalent to a complete star-allele/diplotype call; phenotype assignment from file data requires an explicit validated diplotype/upstream caller and a consistent genotype record.",
        "CYP2D6 copy-number variation, hybrid genes, phasing, and complex structural variation require validated specialized methods.",
        "Risk indicators are drug-specific and evidence-dependent; absence of a matched rule does not imply absence of toxicological risk.",
        "Local laboratory validation, reference materials, SOPs, quality management, privacy controls, ethics approval, and professional review are required before clinical deployment.",
    ]
    for item in limitations:
        story.append(Paragraph(f"• {item}", styles["Small"]))

    story.append(Spacer(1, 10))
    story.append(Paragraph(
        "Evidence framework: FDA pharmacogenetic association table, CPIC guideline resources, and PharmVar allele nomenclature resources. "
        "Evidence sources should be re-reviewed and versioned before each production release.",
        styles["Small"]
    ))

    doc.build(story)
    return buffer.getvalue()


# ============================================================
# SIDEBAR
# ============================================================
with st.sidebar:
    st.header("⚙️ Case Setup")
    case_id = st.text_input("Case ID", value=f"CASE-{datetime.now().strftime('%Y%m%d-%H%M%S')}")
    case_type = st.selectbox(
        "Case type",
        ["Toxicology / overdose assessment", "Addiction-treatment support",
         "Forensic case", "Research / validation case"]
    )

    st.divider()
    st.header("💊 Exposure Context")
    drug = st.selectbox(
        "Drug / substance assessed",
        ["Codeine", "Tramadol", "Clopidogrel", "Amitriptyline",
         "Omeprazole", "Methadone", "Other / not in evidence map"]
    )
    overdose_history = st.selectbox("Previous overdose history", ["Unknown", "No", "Yes"])
    polypharmacy = st.selectbox("Polypharmacy / interacting medicines", ["Unknown", "No", "Yes"])
    acute_exposure = st.selectbox("Current acute exposure", ["Unknown", "No", "Yes"])

    st.divider()
    st.header("🧬 Genomic Input")
    input_method = st.radio(
        "Input method",
        ["VCF / genomic file", "CSV / TXT table", "Interactive panel"]
    )

    st.caption(f"Software v{APP_VERSION} | Build {BUILD_DATE}")


# ============================================================
# MAIN INPUT
# ============================================================
st.markdown('<div class="section-title">1. Genomic Input & QC</div>', unsafe_allow_html=True)

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
        variant_df, parse_warnings = parse_vcf(raw)
    else:
        variant_df, parse_warnings = parse_csv_txt(uploaded)

    qc = qc_variant_table(variant_df)
else:
    qc = {
        "status": "NOT RUN",
        "variant_count": 0,
        "recognized_variants": 0,
        "positive_mapped_variants": 0,
        "missing_genotypes": 0,
        "low_quality": 0,
        "warnings": ["No genomic file supplied. Use the interactive panel or upload validated genomic data."]
    }

qc_cols = st.columns(6)
qc_cols[0].metric("Structural QC", qc["status"])
qc_cols[1].metric("Variants parsed", qc["variant_count"])
qc_cols[2].metric("Mapped loci", qc["recognized_variants"])
qc_cols[3].metric("Positive mapped", qc.get("positive_mapped_variants", 0))
qc_cols[4].metric("Missing GT", qc["missing_genotypes"])
qc_cols[5].metric("Invalid GT", qc.get("invalid_genotypes", 0))
st.caption(f"Interpretation gate: {qc.get('interpretation_status', 'REVIEW REQUIRED')}. Structural QC PASS does not mean analytical or clinical validation.")

if qc["warnings"]:
    with st.expander("QC warnings / review flags", expanded=True):
        for w in qc["warnings"]:
            st.warning(w)

if not variant_df.empty:
    with st.expander("Parsed genomic records", expanded=False):
        st.dataframe(variant_df, use_container_width=True, hide_index=True)

if parse_warnings:
    with st.expander("Parser messages"):
        for w in parse_warnings:
            st.info(w)


# ============================================================
# INTERACTIVE PANEL
# ============================================================
st.markdown('<div class="section-title">2. Phenotype / Genotype Interpretation</div>', unsafe_allow_html=True)

manual_values = {}

if input_method == "Interactive panel":
    st.info(
        "Interactive mode is intended for validated genotype/diplotype entry or research testing. "
        "It does not invent a genotype when data are missing."
    )

    c1, c2 = st.columns(2)
    with c1:
        manual_values["CYP2D6"] = st.selectbox(
            "CYP2D6 diplotype",
            ["Not assessed", "*1/*1", "*1/*2", "*1/*4", "*4/*4", "*1xN"],
        )
        manual_values["CYP2C19"] = st.selectbox(
            "CYP2C19 diplotype",
            ["Not assessed", "*1/*1", "*1/*2", "*2/*2", "*1/*17", "*2/*17"],
        )
        manual_values["CYP3A4"] = st.selectbox(
            "CYP3A4",
            ["Not assessed", "*1/*1", "*1/*22"],
        )
        manual_values["CYP2B6"] = st.selectbox(
            "CYP2B6",
            ["Not assessed", "*1/*1", "*1/*6", "*6/*6"],
        )
        manual_values["CYP2C9"] = st.selectbox(
            "CYP2C9 diplotype",
            ["Not assessed", "*1/*1", "*1/*2", "*1/*3", "*2/*2", "*2/*3", "*3/*3"],
        )

    with c2:
        manual_values["OPRM1"] = st.selectbox(
            "OPRM1",
            ["Not assessed", "A118G variant reported"],
        )
        manual_values["ABCB1"] = st.selectbox(
            "ABCB1",
            ["Not assessed", "rs2032582 variant reported", "rs1045642 variant reported", "Both variants reported"],
        )
        manual_values["HTR2A"] = st.selectbox("HTR2A", ["Not assessed", "rs6311 reported", "rs6313 reported", "Both variants reported"])
        manual_values["ANKK1"] = st.selectbox("ANKK1 / DRD2-associated marker", ["Not assessed", "rs1800497 reported"])
        manual_values["COMT"] = st.selectbox("COMT", ["Not assessed", "rs4680 Val158Met reported"])

    phenotypes, findings = manual_panel_interpretation(manual_values)

else:
    findings = variant_level_findings(variant_df)
    phenotypes = {}

    # Variant-level parsing does not pretend to infer a full diplotype.
    # If a file contains explicit diplotype columns, honor them.
    if not variant_df.empty:
        for gene in ["CYP2D6", "CYP2C19", "CYP3A4", "CYP2B6"]:
            sub = variant_df[variant_df["Gene"].astype(str).str.upper() == gene]
            if sub.empty:
                continue

            # Optional explicit diplotype column from CSV.
            # The supplied diplotype has priority over rsID-only interpretation,
            # but only when QC found no genotype/diplotype conflict for this gene.
            if "Diplotype" in sub.columns:
                values = [normalize_diplotype(v) for v in sub["Diplotype"].tolist() if normalize_diplotype(v)]
                unique_values = list(dict.fromkeys(values))
                gene_conflicts = [c for c in qc.get("diplotype_conflicts", []) if c.get("gene") == gene]
                if gene_conflicts:
                    st.error(f"{gene}: reported diplotype conflicts with genotype data. Phenotype interpretation is blocked for this gene until the source data are resolved.")
                elif len(unique_values) > 1:
                    st.error(f"{gene}: multiple different diplotypes were supplied. Phenotype interpretation is blocked until the source data are resolved.")
                elif unique_values:
                    val = unique_values[0]
                    pheno, reason = phenotype_from_diplotype(gene, val)
                    phenotypes[gene] = (pheno, f"Reported diplotype {val}; {reason}")
                    if pheno == "Unresolved":
                        st.warning(f"{gene}: supplied diplotype {val} is not covered by the conservative built-in phenotype rules. Review is required.")

        if not phenotypes and findings:
            st.info(
                "Variant-level findings were parsed. A full star-allele/diplotype phenotype "
                "is not inferred automatically from isolated variants in this version. "
                "Supply a validated diplotype/caller output for phenotype assignment."
            )

# Optional manual OPRM1 interpretation for file mode is deliberately not defaulted.
manual_oprm1 = "Not assessed"
if input_method != "Interactive panel":
    manual_oprm1 = st.selectbox(
        "Optional OPRM1 interpretation override",
        ["Not assessed", "A118G variant reported"],
        help="Only use when supported by a validated laboratory result."
    )

receptor_status, receptor_reason = receptor_assessment(findings, manual_oprm1)
metabolic_status, metabolic_reason = metabolic_summary(phenotypes)

st.markdown('<div class="section-title">3. Core Results</div>', unsafe_allow_html=True)

r1, r2, r3 = st.columns(3)
r1.metric("Metabolic Status", metabolic_status)
r2.metric("Receptor Assessment", receptor_status)

# Drug-specific risk
exposure_factors = {
    "overdose_history": overdose_history,
    "polypharmacy": polypharmacy,
    "acute_exposure": acute_exposure,
}

if drug == "Other / not in evidence map":
    risk_status = "No validated drug-specific rule available"
    matches = []
    risk_context = "Select a supported drug/substance or add a locally validated evidence rule."
else:
    risk_status, matches, risk_context = risk_engine(
        drug, phenotypes, receptor_status, exposure_factors
    )

r3.metric("Toxicology / PGx Indicator", risk_status)

with st.expander("Explain the core results", expanded=True):
    st.write("**Metabolic interpretation:**", metabolic_reason)
    st.write("**Receptor interpretation:**", receptor_reason)
    st.write("**Risk interpretation:**", risk_context)

# -----------------------------
# Evidence table
# -----------------------------
st.markdown('<div class="section-title">4. Explainable Evidence</div>', unsafe_allow_html=True)

if findings:
    fdf = pd.DataFrame(findings)
    st.dataframe(fdf, use_container_width=True, hide_index=True)
else:
    st.warning(
        "No supported genetic finding has been established. This is not a negative genetic result; "
        "it means the current evidence map cannot establish an interpretable finding."
    )

if matches:
    edf = pd.DataFrame(matches)
    st.dataframe(
        edf[["gene", "indicator", "severity", "mechanism", "evidence"]],
        use_container_width=True,
        hide_index=True
    )
else:
    st.info("No drug-specific evidence rule matched the current phenotype.")

with st.expander("Expanded variant panel — evidence-aware interpretation"):
    st.markdown(
        "The expanded panel records ANKK1/DRD2-associated, HTR2A, ABCB1, CYP2C9 and COMT variants. "
        "Research-level markers are reported as findings and are not converted automatically into a diagnosis, "
        "addiction prediction, or universal overdose label. CYP2C9 findings can contribute to drug-specific metabolic interpretation when the applicable evidence rule is available."
    )
    panel_rows = []
    for rsid, ann in VARIANT_MAP.items():
        if ann.get("gene") in {"ANKK1", "HTR2A", "ABCB1", "CYP2C9", "COMT"}:
            panel_rows.append({"Gene": ann.get("gene"), "Variant": rsid, "Interpretation": ann.get("label"), "Evidence level": ann.get("evidence_level", "Not specified"), "Clinical actionability": "Yes" if ann.get("clinical_actionability") else "No / research"})
    st.dataframe(pd.DataFrame(panel_rows), use_container_width=True, hide_index=True)

# -----------------------------
# Field evaluation variables
# -----------------------------
st.markdown('<div class="section-title">5. Field Evaluation / Case Context</div>', unsafe_allow_html=True)

f1, f2, f3 = st.columns(3)
family_history = f1.selectbox("First-degree family history relevant to substance use / adverse drug response", ["Unknown", "No", "Yes"])
recurrent_overdose = f2.selectbox("Recurrent overdose/adverse-event history", ["Unknown", "No", "Yes"])
toxicology_confirmation = f3.selectbox("Toxicology confirmation available", ["Unknown", "No", "Yes"])

clinical_notes = st.text_area(
    "Professional notes / toxicology findings",
    placeholder="Enter only information permitted by the site's privacy and research SOP."
)

# -----------------------------
# Report generation
# -----------------------------
st.markdown('<div class="section-title">6. Professional Report</div>', unsafe_allow_html=True)

report = {
    "case_id": case_id,
    "case_type": case_type,
    "drug": drug,
    "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "file_hash": file_hash,
    "qc": qc,
    "findings": findings,
    "phenotypes": phenotypes,
    "metabolic_status": metabolic_status,
    "receptor_status": receptor_status,
    "risk_status": risk_status,
    "matches": matches,
    "context": (
        f"Previous overdose history: {overdose_history}; "
        f"polypharmacy: {polypharmacy}; acute exposure: {acute_exposure}; "
        f"first-degree family history: {family_history}; "
        f"recurrent overdose/adverse-event history: {recurrent_overdose}; "
        f"toxicology confirmation: {toxicology_confirmation}. "
        f"Clinical notes: {clinical_notes or 'None provided'}"
    ),
}

if st.button("🧾 Generate PDF Report", type="primary", use_container_width=True):
    if not REPORTLAB_OK:
        st.error("PDF generation requires reportlab. Install it from requirements.txt.")
    else:
        try:
            pdf = make_pdf(report)
            st.success("Report generated. Review the report professionally before any clinical action.")
            st.download_button(
                "⬇️ Download PDF",
                data=pdf,
                file_name=f"{case_id}_PGx_Toxicology_Report.pdf",
                mime="application/pdf",
                use_container_width=True,
            )
        except Exception as exc:
            st.error(f"Could not generate report: {exc}")

# -----------------------------
# Validation / governance panel
# -----------------------------
st.markdown('<div class="section-title">7. Validation & Governance</div>', unsafe_allow_html=True)

with st.expander("What must be validated before hospital/center deployment"):
    st.markdown("""
**Analytical validation**
- Accuracy of VCF/CSV parsing against reference files.
- Genotype concordance against a validated laboratory method.
- QC thresholds and handling of missing/low-quality calls.
- Reproducibility across files/operators/runs.
- CYP2D6 copy-number and structural-variant handling if that gene is clinically reported.

**Clinical / interpretation validation**
- Gene–drug rules must be versioned and linked to an authoritative source.
- Phenotype assignment must be checked against a validated PGx caller/reference dataset.
- Risk indicators must be evaluated against defined clinical/toxicological endpoints.
- False-positive and false-negative behavior must be measured.

**Operational validation**
- User authentication and role control.
- Audit trail.
- Data retention policy.
- De-identification / pseudonymization.
- Backup and incident procedures.
- Local SOP and professional sign-off.

**Important:** This prototype intentionally refuses to call missing genomic data "high risk" and does not treat OPRM1 A118G as a universal overdose predictor.
""")

st.caption(
    "Evidence resources used by the built-in rule set: FDA Pharmacogenetic Associations, "
    "CPIC guideline resources, and PharmVar allele resources. Evidence should be re-checked "
    "and versioned before production release."
)
