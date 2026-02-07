"""
Longevity Genetics v5.0

Complete coverage of:
- FOXO3 (strongest longevity marker)
- APOE (ε2/ε3/ε4 effects)
- TERT/telomere
- CETP (cholesterol)
- IL6 (inflammation)
- KLOTHO (aging)
- SIRT genes
- Centenarian-enriched variants

All markers with PMID references.

NOTE: Longevity is highly polygenic and environmentally influenced.
Genetics explain ~25-30% of lifespan variation.
"""

from enum import Enum
from typing import Dict, List, Any, Optional

class LongevityPotential(Enum):
    HIGH = "high"              # Multiple favorable variants
    ABOVE_AVERAGE = "above_average"
    AVERAGE = "average"
    BELOW_AVERAGE = "below_average"

# =============================================================================
# FOXO3 - STRONGEST LONGEVITY GENE
# =============================================================================

FOXO3_MARKERS = {
    "rs2802292": {
        "gene": "FOXO3",
        "variant": "FOXO3 longevity variant",
        "function": "Forkhead transcription factor - stress resistance, autophagy",
        "risk_allele": "G",  # G is longevity-associated
        "frequency": {"EUR": 0.35, "EAS": 0.30, "AFR": 0.30},
        "effect": {
            "GG": "Longevity-associated genotype - enriched in centenarians",
            "GT": "One longevity allele",
            "TT": "Common genotype"
        },
        "category": "longevity",
        "evidence": "definitive",
        "pmid": ["18765803", "19706745", "25855726"],
        "longevity_effect": {
            "GG": 1.5,  # OR for longevity ~1.5
            "GT": 1.2,
            "TT": 1.0
        },
        "actionable": {
            "GG": [
                "FOXO3 longevity genotype (replicated in multiple populations)",
                "Enriched 2-3x in centenarians vs controls",
                "Associated with healthy aging and stress resistance",
                "Lifestyle factors still most important",
                "This is one of the most replicated longevity variants"
            ]
        }
    },
    "rs2764264": {
        "gene": "FOXO3",
        "variant": "FOXO3 secondary variant",
        "function": "FOXO3 regulation",
        "risk_allele": "C",
        "frequency": {"EUR": 0.35, "EAS": 0.30},
        "effect": "Associated with longevity (LD with rs2802292)",
        "category": "longevity",
        "evidence": "strong",
        "pmid": ["18765803"]
    },
    "rs13217795": {
        "gene": "FOXO3",
        "variant": "FOXO3 variant",
        "function": "FOXO3 activity",
        "risk_allele": "C",
        "frequency": {"EUR": 0.35},
        "effect": "Longevity-associated",
        "category": "longevity",
        "evidence": "moderate",
        "pmid": ["18765803"]
    },
}

# =============================================================================
# APOE - DUAL LONGEVITY AND DISEASE GENE
# =============================================================================

APOE_LONGEVITY = {
    "rs429358": {
        "gene": "APOE",
        "variant": "APOE ε4 determinant",
        "function": "Lipid transport, brain health",
        "risk_allele": "C",  # C = ε4
        "frequency": {"EUR": 0.15, "AFR": 0.25, "EAS": 0.08},
        "category": "longevity",
        "evidence": "definitive",
        "pmid": ["8346443", "19734902"],
        "longevity_effect": {
            "note": "ε4 = reduced longevity, ε2 = increased longevity",
            "e4/e4": 0.7,  # Reduced odds of reaching old age
            "e3/e4": 0.85,
            "e3/e3": 1.0,  # Reference
            "e2/e3": 1.15,  # Increased longevity
            "e2/e2": 1.2
        }
    },
    "rs7412": {
        "gene": "APOE",
        "variant": "APOE ε2 determinant",
        "function": "Combined with rs429358 for APOE genotype",
        "risk_allele": "T",  # T = ε2 (protective)
        "frequency": {"EUR": 0.08},
        "category": "longevity",
        "evidence": "definitive",
        "pmid": ["8346443"],
        "actionable": {
            "e2_carrier": [
                "APOE ε2 is PROTECTIVE for longevity",
                "Lower cardiovascular and Alzheimer's risk",
                "May tolerate saturated fat better"
            ],
            "e4_carrier": [
                "APOE ε4 associated with reduced longevity",
                "Higher Alzheimer's and cardiovascular risk",
                "Lifestyle factors become more important"
            ]
        }
    },
}

# =============================================================================
# TERT/TELOMERE
# =============================================================================

TERT_MARKERS = {
    "rs2736100": {
        "gene": "TERT",
        "variant": "Telomerase reverse transcriptase",
        "function": "Telomere maintenance",
        "risk_allele": "C",
        "frequency": {"EUR": 0.50, "EAS": 0.45},
        "effect": {
            "CC": "Longer telomeres",
            "AA": "Shorter telomeres"
        },
        "category": "longevity",
        "evidence": "strong",
        "pmid": ["21731173", "23446900"],
        "longevity_effect": {
            "CC": 1.1,
            "CA": 1.0,
            "AA": 0.95
        },
        "actionable": {
            "CC": [
                "Associated with longer telomeres",
                "Telomere length associated with biological aging",
                "Still highly influenced by lifestyle (exercise, stress, diet)"
            ]
        }
    },
    "rs7726159": {
        "gene": "TERT",
        "variant": "TERT promoter variant",
        "function": "TERT expression",
        "risk_allele": "C",
        "frequency": {"EUR": 0.30},
        "effect": "Telomere length association",
        "category": "longevity",
        "evidence": "moderate",
        "pmid": ["23446900"]
    },
}

# =============================================================================
# CETP - CHOLESTEROL AND LONGEVITY
# =============================================================================

CETP_MARKERS = {
    "rs5882": {
        "gene": "CETP",
        "variant": "I405V",
        "function": "Cholesteryl ester transfer protein",
        "risk_allele": "A",  # Val = lower CETP
        "frequency": {"EUR": 0.35, "EAS": 0.50},
        "effect": {
            "AA": "Val/Val - Lower CETP activity, higher HDL",
            "note": "Enriched in centenarians"
        },
        "category": "longevity",
        "evidence": "strong",
        "pmid": ["15655135", "17848700"],
        "longevity_effect": {
            "AA": 1.2
        },
        "actionable": {
            "AA": [
                "Lower CETP activity - higher HDL cholesterol",
                "Favorable lipoprotein profile",
                "Enriched in centenarians",
                "Note: CETP inhibitor drugs did not prevent CVD"
            ]
        }
    },
    "rs1800775": {
        "gene": "CETP",
        "variant": "TaqIB (-629C>A)",
        "function": "CETP expression",
        "risk_allele": "A",  # Lower CETP
        "frequency": {"EUR": 0.42},
        "effect": {
            "AA": "B2/B2 - Higher HDL, lower CETP"
        },
        "category": "longevity",
        "evidence": "strong",
        "pmid": ["17848700", "15655135"]
    },
}

# =============================================================================
# IL6 - INFLAMMATION AND AGING
# =============================================================================

IL6_LONGEVITY = {
    "rs1800795": {
        "gene": "IL6",
        "variant": "-174G>C",
        "function": "Interleukin-6 expression",
        "risk_allele": "C",
        "frequency": {"EUR": 0.45},
        "effect": {
            "GG": "Lower IL-6 - less inflammation, better longevity",
            "CC": "Higher IL-6 - more inflammation"
        },
        "category": "longevity",
        "evidence": "strong",
        "pmid": ["12730251", "18651211"],
        "longevity_effect": {
            "GG": 1.15,
            "GC": 1.0,
            "CC": 0.9
        },
        "actionable": {
            "GG": [
                "Lower baseline inflammation",
                "Associated with better longevity outcomes",
                "Anti-inflammatory lifestyle amplifies benefit"
            ],
            "CC": [
                "Higher baseline inflammation",
                "Anti-inflammatory diet/lifestyle especially important",
                "Omega-3, exercise, stress management beneficial"
            ]
        }
    },
}

# =============================================================================
# KLOTHO - AGING SUPPRESSOR
# =============================================================================

KLOTHO_MARKERS = {
    "rs9536314": {
        "gene": "KL",
        "variant": "F352V (Klotho-VS)",
        "function": "Klotho aging suppressor protein",
        "risk_allele": "T",  # VS variant
        "frequency": {"EUR": 0.15, "AFR": 0.10, "EAS": 0.05},
        "effect": {
            "GT": "KL-VS heterozygote - associated with longevity and cognition",
            "note": "Heterozygotes advantaged; homozygotes may have reduced effect"
        },
        "category": "longevity",
        "evidence": "strong",
        "pmid": ["17570871", "24630764"],
        "longevity_effect": {
            "GT": 1.15,
            "GG": 1.0,
            "TT": 0.95  # Homozygotes may not have same benefit
        },
        "actionable": {
            "GT": [
                "Klotho-VS heterozygote - favorable for aging",
                "Associated with better cognitive aging",
                "Higher serum Klotho levels"
            ]
        }
    },
    "rs9527025": {
        "gene": "KL",
        "variant": "Klotho C370S",
        "function": "In cis with F352V",
        "risk_allele": "G",
        "frequency": {"EUR": 0.15},
        "effect": "Part of KL-VS haplotype",
        "category": "longevity",
        "evidence": "strong",
        "pmid": ["17570871"]
    },
}

# =============================================================================
# SIRT GENES - SIRTUINS
# =============================================================================

SIRT_MARKERS = {
    "rs3758391": {
        "gene": "SIRT1",
        "variant": "SIRT1 promoter",
        "function": "Sirtuin 1 - NAD+ dependent deacetylase",
        "risk_allele": "C",
        "frequency": {"EUR": 0.35, "EAS": 0.25},
        "effect": {
            "CC": "Associated with longevity in some studies"
        },
        "category": "longevity",
        "evidence": "moderate",
        "pmid": ["18048405", "21106703"],
        "longevity_effect": {
            "CC": 1.1
        },
        "actionable": {
            "CC": [
                "May have enhanced SIRT1 function",
                "Caloric restriction and exercise activate SIRT1",
                "Resveratrol claimed to activate (evidence mixed)"
            ]
        }
    },
    "rs10823108": {
        "gene": "SIRT3",
        "variant": "SIRT3 enhancer",
        "function": "Mitochondrial sirtuin",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30},
        "effect": "Associated with longevity",
        "category": "longevity",
        "evidence": "moderate",
        "pmid": ["21211615"],
        "longevity_effect": {
            "AA": 1.1
        }
    },
}

# =============================================================================
# ADDITIONAL LONGEVITY-ASSOCIATED MARKERS
# =============================================================================

OTHER_LONGEVITY_MARKERS = {
    "rs1800795": {
        "gene": "IL6",
        "reference": "See IL6_LONGEVITY"
    },
    "rs2075650": {
        "gene": "TOMM40/APOE region",
        "variant": "TOMM40 variant",
        "function": "Mitochondrial protein import (near APOE)",
        "risk_allele": "G",
        "frequency": {"EUR": 0.15},
        "effect": "Associated with APOE effects and longevity",
        "category": "longevity",
        "evidence": "moderate",
        "pmid": ["20708966", "23357107"],
        "note": "In strong LD with APOE - may not be independent"
    },
    "rs1801131": {
        "gene": "MTHFR",
        "variant": "A1298C",
        "function": "Methylation - DNA repair",
        "risk_allele": "C",
        "frequency": {"EUR": 0.33},
        "effect": "Some longevity studies show association",
        "category": "longevity",
        "evidence": "moderate",
        "pmid": ["15059627"]
    },
    "rs4402960": {
        "gene": "IGF2BP2",
        "variant": "IGF pathway variant",
        "function": "Insulin/IGF signaling",
        "risk_allele": "T",
        "frequency": {"EUR": 0.30},
        "effect": {
            "note": "Lower IGF signaling associated with longevity in model organisms"
        },
        "category": "longevity",
        "evidence": "moderate",
        "pmid": ["17463246"]
    },
}

# =============================================================================
# COMBINE ALL LONGEVITY MARKERS
# =============================================================================

LONGEVITY_COMPREHENSIVE_MARKERS = {
    **FOXO3_MARKERS,
    **APOE_LONGEVITY,
    **TERT_MARKERS,
    **CETP_MARKERS,
    **IL6_LONGEVITY,
    **KLOTHO_MARKERS,
    **SIRT_MARKERS,
    **OTHER_LONGEVITY_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def calculate_longevity_score(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate genetic longevity potential score."""
    
    score = 1.0  # Multiplicative odds ratios
    factors = []
    
    # FOXO3 (strongest)
    foxo3 = genotypes.get("rs2802292", "TT")
    if foxo3 == "GG":
        score *= 1.5
        factors.append("FOXO3 GG - strongest longevity variant (OR 1.5)")
    elif foxo3 == "GT":
        score *= 1.2
        factors.append("FOXO3 GT - one longevity allele")
    
    # APOE
    # Simplified - would need full APOE typing
    rs429358 = genotypes.get("rs429358", "TT")
    rs7412 = genotypes.get("rs7412", "CC")
    
    if rs429358 == "CC":  # ε4/ε4 or ε3/ε4
        score *= 0.75
        factors.append("APOE ε4 carrier - reduced longevity association")
    elif rs7412 == "TT":  # ε2/ε2
        score *= 1.15
        factors.append("APOE ε2 homozygous - longevity protective")
    elif rs7412 == "CT":  # ε2/ε3
        score *= 1.1
        factors.append("APOE ε2 carrier - longevity protective")
    
    # CETP
    cetp = genotypes.get("rs5882", "GG")
    if cetp == "AA":
        score *= 1.1
        factors.append("CETP I405V Val/Val - higher HDL, centenarian-enriched")
    
    # IL6
    il6 = genotypes.get("rs1800795", "GC")
    if il6 == "GG":
        score *= 1.1
        factors.append("IL6 GG - lower inflammation")
    elif il6 == "CC":
        score *= 0.95
        factors.append("IL6 CC - higher inflammation")
    
    # KLOTHO
    klotho = genotypes.get("rs9536314", "GG")
    if klotho == "GT":
        score *= 1.1
        factors.append("Klotho-VS heterozygote - favorable aging")
    
    # TERT
    tert = genotypes.get("rs2736100", "CA")
    if tert == "CC":
        score *= 1.05
        factors.append("TERT CC - longer telomeres")
    
    # Determine potential
    if score >= 1.4:
        potential = LongevityPotential.HIGH
    elif score >= 1.15:
        potential = LongevityPotential.ABOVE_AVERAGE
    elif score >= 0.9:
        potential = LongevityPotential.AVERAGE
    else:
        potential = LongevityPotential.BELOW_AVERAGE
    
    return {
        "longevity_score": round(score, 2),
        "potential": potential.value,
        "favorable_factors": factors,
        "interpretation": f"Genetic longevity score: {round(score, 2)}x population average",
        "important_note": "Genetics explain only ~25-30% of lifespan. Lifestyle factors dominate.",
        "lifestyle_recommendations": [
            "Regular exercise (most impactful modifiable factor)",
            "Mediterranean or whole-food diet",
            "Don't smoke",
            "Moderate alcohol or none",
            "Maintain healthy weight",
            "Social connections and purpose",
            "Quality sleep",
            "Manage chronic stress"
        ],
        "pmid": ["18765803", "8346443", "15655135"]
    }

def generate_longevity_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate comprehensive longevity genetics report."""
    
    longevity_score = calculate_longevity_score(genotypes)
    
    # Key markers summary
    foxo3 = genotypes.get("rs2802292", "Unknown")
    
    return {
        **longevity_score,
        "key_markers": {
            "FOXO3": foxo3,
            "APOE_ε4": "carrier" if genotypes.get("rs429358") == "CC" else "non-carrier",
            "CETP_favorable": genotypes.get("rs5882") == "AA",
            "IL6_low_inflammation": genotypes.get("rs1800795") == "GG",
        },
        "markers_analyzed": sum(1 for rs in LONGEVITY_COMPREHENSIVE_MARKERS if rs in genotypes),
    }

# Export
__all__ = [
    'LONGEVITY_COMPREHENSIVE_MARKERS',
    'FOXO3_MARKERS',
    'APOE_LONGEVITY',
    'TERT_MARKERS',
    'CETP_MARKERS',
    'IL6_LONGEVITY',
    'KLOTHO_MARKERS',
    'SIRT_MARKERS',
    'OTHER_LONGEVITY_MARKERS',
    'LongevityPotential',
    'calculate_longevity_score',
    'generate_longevity_report',
]
