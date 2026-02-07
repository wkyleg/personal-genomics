"""
Medical Special Categories v5.0

Critical medical genetics:
- Anesthesia (BCHE, RYR1 for malignant hyperthermia)
- Blood type inference (ABO, RH)
- Infection resistance (CCR5-delta32 HIV, FUT2 norovirus, DARC malaria)
- Autoimmune HLA panel
- Celiac certainty (DQ2.5, DQ8, DQ2.2 rule-out)
- Vaccine response (HLA associations)

All markers with PMID references and clinical implications.
"""

from enum import Enum
from typing import Dict, List, Any, Optional

class AnesthesiaRisk(Enum):
    CRITICAL = "critical"      # Life-threatening risk
    ELEVATED = "elevated"      # Needs precautions
    NORMAL = "normal"

class BloodType(Enum):
    O_POS = "O+"
    O_NEG = "O-"
    A_POS = "A+"
    A_NEG = "A-"
    B_POS = "B+"
    B_NEG = "B-"
    AB_POS = "AB+"
    AB_NEG = "AB-"
    UNKNOWN = "unknown"

# =============================================================================
# ANESTHESIA - MALIGNANT HYPERTHERMIA
# =============================================================================

ANESTHESIA_MH_MARKERS = {
    "rs121918592": {
        "gene": "RYR1",
        "variant": "R614C",
        "function": "Ryanodine receptor - calcium release in muscle",
        "risk_allele": "T",
        "frequency": {"EUR": 0.0001},  # Rare but critical
        "effect": {
            "CT/TT": "Malignant hyperthermia susceptibility (MHS)"
        },
        "category": "anesthesia",
        "condition": "Malignant Hyperthermia",
        "severity": AnesthesiaRisk.CRITICAL,
        "evidence": "definitive",
        "pmid": ["21880859", "9762081"],
        "actionable": {
            "carrier": [
                "⚠️ MALIGNANT HYPERTHERMIA SUSCEPTIBLE",
                "AVOID succinylcholine and ALL volatile anesthetics",
                "Use total intravenous anesthesia (TIVA) only",
                "Dantrolene must be immediately available",
                "Medical alert bracelet REQUIRED",
                "Inform ALL healthcare providers",
                "Family members should be tested",
                "70% mortality if untreated"
            ]
        }
    },
    "rs28933396": {
        "gene": "RYR1",
        "variant": "G341R",
        "function": "RYR1 MH variant",
        "risk_allele": "A",
        "frequency": {"EUR": 0.0001},
        "effect": "Malignant hyperthermia susceptibility",
        "category": "anesthesia",
        "severity": AnesthesiaRisk.CRITICAL,
        "evidence": "definitive",
        "pmid": ["21880859"]
    },
    "rs121918593": {
        "gene": "RYR1",
        "variant": "G2434R",
        "function": "RYR1 MH variant",
        "risk_allele": "A",
        "frequency": {"EUR": 0.0001},
        "effect": "Malignant hyperthermia susceptibility",
        "category": "anesthesia",
        "severity": AnesthesiaRisk.CRITICAL,
        "evidence": "definitive",
        "pmid": ["21880859"]
    },
    "rs118192178": {
        "gene": "RYR1",
        "variant": "R163C",
        "function": "RYR1 MH variant",
        "risk_allele": "T",
        "frequency": {"EUR": 0.0001},
        "effect": "MH susceptibility",
        "category": "anesthesia",
        "severity": AnesthesiaRisk.CRITICAL,
        "evidence": "definitive",
        "pmid": ["21880859"]
    },
    "rs772226819": {
        "gene": "CACNA1S",
        "variant": "R1086H",
        "function": "Calcium channel - secondary MH gene",
        "risk_allele": "A",
        "frequency": {"EUR": 0.00005},
        "effect": "MH susceptibility (1-2% of MH cases)",
        "category": "anesthesia",
        "severity": AnesthesiaRisk.CRITICAL,
        "evidence": "definitive",
        "pmid": ["21880859"]
    },
}

# =============================================================================
# ANESTHESIA - PSEUDOCHOLINESTERASE DEFICIENCY
# =============================================================================

BCHE_MARKERS = {
    "rs1799807": {
        "gene": "BCHE",
        "variant": "Atypical (D70G)",
        "function": "Butyrylcholinesterase - hydrolyzes succinylcholine",
        "risk_allele": "A",
        "frequency": {"EUR": 0.02, "SAS": 0.04},
        "effect": {
            "AA": "Pseudocholinesterase deficiency - prolonged paralysis (2-8 hours)",
            "AG": "Carrier - may have mild prolongation"
        },
        "category": "anesthesia",
        "severity": AnesthesiaRisk.ELEVATED,
        "evidence": "definitive",
        "pmid": ["1349607", "8630723"],
        "actionable": {
            "AA": [
                "Pseudocholinesterase deficiency",
                "Succinylcholine will cause PROLONGED paralysis (2-8 hours)",
                "Inform anesthesiologist BEFORE any surgery",
                "Use rocuronium + sugammadex as alternative",
                "Not dangerous if recognized - just prolonged ventilation needed",
                "Family members should be tested"
            ],
            "AG": [
                "Carrier - may have slightly prolonged response to succinylcholine",
                "Inform anesthesiologist"
            ]
        }
    },
    "rs28933389": {
        "gene": "BCHE",
        "variant": "K variant",
        "function": "Reduced BCHE activity (~66%)",
        "risk_allele": "A",
        "frequency": {"EUR": 0.12},
        "effect": {
            "AA": "Mild BCHE reduction - may prolong succinylcholine effect if compound het"
        },
        "category": "anesthesia",
        "severity": AnesthesiaRisk.ELEVATED,
        "evidence": "moderate",
        "pmid": ["1349607"]
    },
    "rs1803274": {
        "gene": "BCHE",
        "variant": "K variant 2",
        "function": "In LD with K variant",
        "risk_allele": "A",
        "frequency": {"EUR": 0.12},
        "effect": "BCHE reduction marker",
        "category": "anesthesia",
        "evidence": "moderate",
        "pmid": ["1349607"]
    },
}

# =============================================================================
# BLOOD TYPE INFERENCE (ABO + RH)
# =============================================================================

ABO_MARKERS = {
    "rs8176719": {
        "gene": "ABO",
        "variant": "O allele determinant (261delG)",
        "function": "Frameshift creating O allele",
        "risk_allele": "del",  # Deletion = O
        "frequency": {"EUR": 0.40, "AFR": 0.50, "EAS": 0.35, "AMR": 0.55},
        "effect": {
            "del/del": "Blood type O",
            "ins/del": "A or B carrier of O",
            "ins/ins": "A and/or B"
        },
        "category": "blood_type",
        "evidence": "definitive",
        "pmid": ["2319781", "1737852"]
    },
    "rs8176746": {
        "gene": "ABO",
        "variant": "B allele marker",
        "function": "Distinguishes A from B",
        "risk_allele": "G",  # G = B, C = A
        "frequency": {"EUR": 0.08, "AFR": 0.15, "EAS": 0.20, "SAS": 0.25},
        "effect": {
            "G": "B allele",
            "C": "A allele"
        },
        "category": "blood_type",
        "evidence": "definitive",
        "pmid": ["2319781"]
    },
    "rs8176747": {
        "gene": "ABO",
        "variant": "A/B modifier",
        "function": "Secondary A/B marker",
        "risk_allele": "G",
        "frequency": {"EUR": 0.10},
        "effect": "Helps distinguish A/B subtypes",
        "category": "blood_type",
        "evidence": "strong",
        "pmid": ["2319781"]
    },
}

RH_MARKERS = {
    "rs590787": {
        "gene": "RHD",
        "variant": "RhD proxy SNP",
        "function": "RhD blood group",
        "risk_allele": "A",  # Presence indicator
        "frequency": {"EUR": 0.85, "AFR": 0.95, "EAS": 0.99},
        "effect": {
            "note": "RhD is determined by gene deletion/presence",
            "proxy": "SNPs can approximate but not definitively type Rh"
        },
        "category": "blood_type",
        "evidence": "moderate",
        "pmid": ["1552912"],
        "note": "Rh typing from SNP arrays is approximate - confirm with serology"
    },
}

# =============================================================================
# INFECTION RESISTANCE
# =============================================================================

INFECTION_RESISTANCE_MARKERS = {
    # CCR5 Delta32 - HIV resistance
    "rs333": {
        "gene": "CCR5",
        "variant": "Delta32 (32bp deletion)",
        "function": "HIV co-receptor - required for HIV entry",
        "risk_allele": "del",
        "frequency": {"EUR": 0.10, "AFR": 0.01, "EAS": 0.01},
        "effect": {
            "del/del": "RESISTANT to R5-tropic HIV infection (most HIV strains)",
            "ins/del": "Slower HIV progression if infected, partial protection"
        },
        "category": "infection",
        "condition": "HIV resistance",
        "evidence": "definitive",
        "pmid": ["8791590", "9399903"],
        "actionable": {
            "del/del": [
                "Natural resistance to most HIV strains (CCR5-tropic)",
                "NOT resistant to X4-tropic HIV strains (rare)",
                "Still use protection - other STIs possible",
                "This is the mutation Berlin/London patients had",
                "~1% of Europeans are homozygous"
            ],
            "ins/del": [
                "Partial HIV protection",
                "Slower progression to AIDS if infected",
                "~10% of Europeans are heterozygous"
            ]
        }
    },
    
    # FUT2 - Norovirus resistance
    "rs601338": {
        "gene": "FUT2",
        "variant": "W154X (Non-secretor)",
        "function": "Fucosyltransferase 2 - gut surface glycans",
        "risk_allele": "A",  # A = non-secretor
        "frequency": {"EUR": 0.45, "AFR": 0.35, "EAS": 0.15},
        "effect": {
            "AA": "Non-secretor - RESISTANT to most norovirus strains",
            "GA/GG": "Secretor - susceptible to norovirus"
        },
        "category": "infection",
        "condition": "Norovirus resistance",
        "evidence": "definitive",
        "pmid": ["12692541", "19578716"],
        "actionable": {
            "AA": [
                "Non-secretor - resistant to MOST norovirus strains",
                "May still get some rarer strains",
                "Won't get the 'cruise ship' norovirus",
                "~20% of Caucasians are non-secretors",
                "Also: lower B12 absorption, different gut microbiome"
            ],
            "GG": [
                "Secretor - susceptible to norovirus",
                "No special protection from food poisoning viruses"
            ]
        }
    },
    
    # DARC - Malaria resistance (Duffy antigen)
    "rs2814778": {
        "gene": "ACKR1/DARC",
        "variant": "Duffy null (GATA box)",
        "function": "Duffy antigen receptor for chemokines",
        "risk_allele": "C",  # C = Duffy null
        "frequency": {"EUR": 0.001, "AFR": 0.95, "AFR_WA": 0.99},
        "effect": {
            "CC": "Duffy null - RESISTANT to P. vivax malaria",
            "note": "Nearly fixed in sub-Saharan Africa due to malaria selection"
        },
        "category": "infection",
        "condition": "Plasmodium vivax malaria resistance",
        "evidence": "definitive",
        "pmid": ["8090753", "11564488"],
        "actionable": {
            "CC": [
                "Duffy null - resistant to P. vivax malaria",
                "NOT resistant to P. falciparum (deadly form)",
                "Still need malaria prophylaxis in endemic areas",
                "Nearly universal in West African ancestry"
            ]
        }
    },
    
    # HBB - Sickle cell trait (malaria protection)
    "rs334": {
        "gene": "HBB",
        "variant": "HbS (E6V) Sickle cell",
        "function": "Beta-globin - sickle cell",
        "risk_allele": "T",  # T = HbS
        "frequency": {"EUR": 0.001, "AFR": 0.10, "AFR_WA": 0.20},
        "effect": {
            "AT": "Sickle cell TRAIT - partial P. falciparum malaria protection",
            "TT": "Sickle cell DISEASE - avoid unless discussing with doctor"
        },
        "category": "infection",
        "condition": "Malaria resistance (sickle trait)",
        "evidence": "definitive",
        "pmid": ["11919001"],
        "actionable": {
            "AT": [
                "Sickle cell trait (carrier)",
                "90% protection against severe falciparum malaria",
                "Rarely causes problems but inform doctors",
                "Avoid extreme exertion at altitude",
                "50% chance to pass to children"
            ],
            "TT": [
                "⚠️ SICKLE CELL DISEASE markers detected",
                "Confirm with hemoglobin electrophoresis",
                "Major medical implications if confirmed"
            ]
        }
    },
}

# =============================================================================
# CELIAC HLA PANEL (COMPLETE)
# =============================================================================

CELIAC_HLA_COMPLETE = {
    "rs2187668": {
        "gene": "HLA-DQA1*05:01",
        "variant": "DQ2.5 alpha chain",
        "function": "Part of HLA-DQ2.5 heterodimer",
        "risk_allele": "T",
        "frequency": {"EUR": 0.20, "AFR": 0.10, "EAS": 0.05},
        "category": "celiac",
        "dq_type": "DQ2.5_alpha",
        "evidence": "definitive",
        "pmid": ["18509540", "20190752", "24867074"]
    },
    "rs7454108": {
        "gene": "HLA-DQB1*02",
        "variant": "DQ2.5 beta chain",
        "function": "Part of HLA-DQ2.5 heterodimer",
        "risk_allele": "C",
        "frequency": {"EUR": 0.25},
        "category": "celiac",
        "dq_type": "DQ2.5_beta",
        "evidence": "definitive",
        "pmid": ["18509540"]
    },
    "rs7775228": {
        "gene": "HLA-DQB1*03:02",
        "variant": "DQ8 marker",
        "function": "HLA-DQ8 - secondary celiac haplotype",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15, "AMR": 0.25},
        "category": "celiac",
        "dq_type": "DQ8",
        "evidence": "definitive",
        "pmid": ["18509540", "20190752"]
    },
    "rs2395182": {
        "gene": "HLA-DQA1*02:01",
        "variant": "DQ2.2 alpha component",
        "function": "Part of DQ2.2 haplotype",
        "risk_allele": "T",
        "frequency": {"EUR": 0.25},
        "category": "celiac",
        "dq_type": "DQ2.2_alpha",
        "evidence": "strong",
        "pmid": ["24867074"]
    },
}

# =============================================================================
# VACCINE RESPONSE
# =============================================================================

VACCINE_RESPONSE_MARKERS = {
    # Hepatitis B vaccine non-response
    "rs2856718": {
        "gene": "HLA-DQ region",
        "variant": "HBV vaccine response",
        "function": "HLA class II antigen presentation",
        "risk_allele": "G",
        "frequency": {"EUR": 0.30},
        "effect": {
            "GG": "May have reduced hepatitis B vaccine response"
        },
        "category": "vaccine",
        "evidence": "moderate",
        "pmid": ["18200439", "25436858"],
        "actionable": {
            "GG": [
                "May not respond well to hepatitis B vaccine",
                "Check anti-HBs titers after vaccination",
                "May need additional booster doses"
            ]
        }
    },
    
    # Measles vaccine response
    "rs2844580": {
        "gene": "HLA region",
        "variant": "Measles antibody response",
        "function": "Immune response to measles antigen",
        "risk_allele": "A",
        "frequency": {"EUR": 0.20},
        "effect": "Associated with measles vaccine response variation",
        "category": "vaccine",
        "evidence": "moderate",
        "pmid": ["21320523"]
    },
    
    # Rubella vaccine response
    "rs2269423": {
        "gene": "TRIM22/TRIM5",
        "variant": "Rubella antibody response",
        "function": "Innate immune response",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35},
        "effect": "Associated with rubella vaccine antibody levels",
        "category": "vaccine",
        "evidence": "moderate",
        "pmid": ["21320523"]
    },
    
    # Influenza vaccine response
    "rs2857149": {
        "gene": "HLA-DQB1",
        "variant": "Influenza vaccine response",
        "function": "Antigen presentation",
        "risk_allele": "C",
        "frequency": {"EUR": 0.25},
        "effect": "Associated with influenza vaccine response",
        "category": "vaccine",
        "evidence": "moderate",
        "pmid": ["23684983"]
    },
}

# =============================================================================
# COMBINE ALL MEDICAL SPECIAL MARKERS
# =============================================================================

MEDICAL_SPECIAL_MARKERS = {
    **ANESTHESIA_MH_MARKERS,
    **BCHE_MARKERS,
    **ABO_MARKERS,
    **RH_MARKERS,
    **INFECTION_RESISTANCE_MARKERS,
    **CELIAC_HLA_COMPLETE,
    **VACCINE_RESPONSE_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def check_anesthesia_alerts(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Check for critical anesthesia-related variants."""
    
    alerts = []
    
    # Check RYR1 for MH
    ryr1_variants = ["rs121918592", "rs28933396", "rs121918593", "rs118192178"]
    for rs in ryr1_variants:
        geno = genotypes.get(rs, "")
        if geno and any(allele in geno for allele in ["T", "A"]):
            alerts.append({
                "gene": "RYR1",
                "severity": "CRITICAL",
                "condition": "Malignant Hyperthermia Susceptibility",
                "message": "AVOID succinylcholine and volatile anesthetics",
                "action": "Medical alert bracelet required. Inform all healthcare providers.",
                "pmid": "21880859"
            })
            break
    
    # Check CACNA1S
    cacna1s = genotypes.get("rs772226819", "GG")
    if "A" in cacna1s:
        alerts.append({
            "gene": "CACNA1S",
            "severity": "CRITICAL",
            "condition": "Malignant Hyperthermia Susceptibility",
            "message": "Secondary MH gene - same precautions as RYR1",
            "pmid": "21880859"
        })
    
    # Check BCHE
    bche = genotypes.get("rs1799807", "GG")
    if bche == "AA":
        alerts.append({
            "gene": "BCHE",
            "severity": "ELEVATED",
            "condition": "Pseudocholinesterase Deficiency",
            "message": "Succinylcholine will cause prolonged paralysis (2-8 hours)",
            "action": "Use rocuronium + sugammadex instead",
            "pmid": "1349607"
        })
    elif bche == "AG":
        alerts.append({
            "gene": "BCHE",
            "severity": "MODERATE",
            "condition": "BCHE Carrier",
            "message": "May have slightly prolonged response to succinylcholine",
            "pmid": "1349607"
        })
    
    return {
        "alerts": alerts,
        "has_critical": any(a["severity"] == "CRITICAL" for a in alerts),
        "recommendation": "Share with anesthesiologist before any procedure" if alerts else "No anesthesia alerts detected"
    }

def infer_blood_type(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Infer ABO blood type from genetics (approximate)."""
    
    # This is simplified - true ABO typing is complex
    rs8176719 = genotypes.get("rs8176719", "ins")
    rs8176746 = genotypes.get("rs8176746", "C")
    
    # Determine ABO
    if "del" in rs8176719 or rs8176719 == "GG":  # O indicator varies by array
        # Check if homozygous O
        if rs8176719 in ["del/del", "GG"]:
            abo = "O"
        else:
            # Heterozygous - check A vs B
            if rs8176746 == "G":
                abo = "B"
            else:
                abo = "A"
    else:
        # No O allele
        if rs8176746 == "GG":
            abo = "B"
        elif rs8176746 == "CC":
            abo = "A"
        elif rs8176746 == "GC":
            abo = "AB"
        else:
            abo = "A"  # Default
    
    return {
        "inferred_abo": abo,
        "confidence": "moderate",
        "note": "SNP-based inference is approximate. Confirm with serology for medical use.",
        "pmid": ["2319781"]
    }

def check_celiac_hla(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Determine celiac HLA risk status."""
    
    dq2_5_alpha = genotypes.get("rs2187668", "CC")
    dq2_5_beta = genotypes.get("rs7454108", "TT")
    dq8 = genotypes.get("rs7775228", "TT")
    dq2_2_alpha = genotypes.get("rs2395182", "CC")
    
    # DQ2.5 (highest risk - 95% of celiac)
    has_dq2_5 = dq2_5_alpha == "TT" or (dq2_5_alpha == "CT" and dq2_5_beta in ["CC", "CT"])
    
    # DQ8 (5-10% of celiac without DQ2)
    has_dq8 = dq8 in ["CC", "CT"]
    
    # DQ2.2 (lower risk)
    has_dq2_2 = dq2_2_alpha in ["TT", "CT"]
    
    # Determine risk
    if has_dq2_5:
        if dq2_5_alpha == "TT":
            risk = "VERY HIGH"
            interpretation = "HLA-DQ2.5 homozygous - highest celiac risk (~10% lifetime)"
        else:
            risk = "HIGH"
            interpretation = "HLA-DQ2.5 heterozygous - elevated celiac risk (~3-5%)"
    elif has_dq8:
        risk = "MODERATE"
        interpretation = "HLA-DQ8 positive - celiac possible (~1-2% risk)"
    elif has_dq2_2:
        risk = "LOW"
        interpretation = "HLA-DQ2.2 only - low but not zero celiac risk"
    else:
        risk = "VERY LOW"
        interpretation = "No DQ2 or DQ8 - celiac essentially ruled out (<0.1%)"
    
    can_rule_out = risk == "VERY LOW"
    
    return {
        "risk_level": risk,
        "interpretation": interpretation,
        "hla_dq2_5": has_dq2_5,
        "hla_dq8": has_dq8,
        "hla_dq2_2": has_dq2_2,
        "celiac_can_be_ruled_out": can_rule_out,
        "recommendation": "tTG-IgA testing if symptomatic" if not can_rule_out else "Celiac very unlikely based on genetics",
        "note": "Do NOT start gluten-free diet before testing",
        "pmid": ["18509540", "20190752"]
    }

def check_infection_resistance(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Check infection resistance variants."""
    
    resistances = []
    
    # CCR5-delta32 (HIV)
    ccr5 = genotypes.get("rs333", "ins/ins")
    if "del" in ccr5:
        if ccr5 == "del/del":
            resistances.append({
                "pathogen": "HIV (R5-tropic strains)",
                "status": "RESISTANT",
                "note": "Homozygous CCR5-delta32 - natural HIV resistance"
            })
        else:
            resistances.append({
                "pathogen": "HIV",
                "status": "PARTIAL PROTECTION",
                "note": "Heterozygous CCR5-delta32 - slower progression"
            })
    
    # FUT2 (Norovirus)
    fut2 = genotypes.get("rs601338", "GG")
    if fut2 == "AA":
        resistances.append({
            "pathogen": "Norovirus (most strains)",
            "status": "RESISTANT",
            "note": "Non-secretor - resistant to common norovirus"
        })
    
    # DARC (P. vivax malaria)
    darc = genotypes.get("rs2814778", "TT")
    if darc == "CC":
        resistances.append({
            "pathogen": "Plasmodium vivax malaria",
            "status": "RESISTANT",
            "note": "Duffy null - common in African ancestry"
        })
    
    # Sickle cell trait (P. falciparum)
    hbs = genotypes.get("rs334", "AA")
    if hbs == "AT":
        resistances.append({
            "pathogen": "Plasmodium falciparum malaria",
            "status": "PARTIAL PROTECTION",
            "note": "Sickle cell trait provides ~90% protection from severe malaria"
        })
    
    return {
        "resistances": resistances,
        "count": len(resistances),
        "pmid": ["8791590", "12692541", "8090753", "11919001"]
    }

def generate_medical_special_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate complete medical special report."""
    
    return {
        "anesthesia": check_anesthesia_alerts(genotypes),
        "blood_type": infer_blood_type(genotypes),
        "celiac_hla": check_celiac_hla(genotypes),
        "infection_resistance": check_infection_resistance(genotypes),
        "markers_analyzed": sum(1 for rs in MEDICAL_SPECIAL_MARKERS if rs in genotypes),
        "critical_disclaimer": "Confirm critical findings with clinical testing. This is screening, not diagnosis."
    }

# Export
__all__ = [
    'MEDICAL_SPECIAL_MARKERS',
    'ANESTHESIA_MH_MARKERS',
    'BCHE_MARKERS',
    'ABO_MARKERS',
    'RH_MARKERS',
    'INFECTION_RESISTANCE_MARKERS',
    'CELIAC_HLA_COMPLETE',
    'VACCINE_RESPONSE_MARKERS',
    'AnesthesiaRisk',
    'BloodType',
    'check_anesthesia_alerts',
    'infer_blood_type',
    'check_celiac_hla',
    'check_infection_resistance',
    'generate_medical_special_report',
]
