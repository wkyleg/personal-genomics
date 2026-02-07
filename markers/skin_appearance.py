"""
Skin & Appearance Genetics v5.0

Complete coverage of:
- UV sensitivity (full MC1R, IRF4, TYR, OCA2, SLC45A2)
- Sunburn and tanning
- Melanoma risk
- Freckling
- Skin aging (MMP1, COL1A1)
- Acne, rosacea, psoriasis, eczema
- Hair color, loss, gray timing, texture
- Eye color prediction

All markers with PMID references.
"""

from enum import Enum
from typing import Dict, List, Any, Optional

class SkinType(Enum):
    TYPE_I = "type_i"     # Always burns, never tans (very fair)
    TYPE_II = "type_ii"   # Usually burns, tans minimally
    TYPE_III = "type_iii" # Sometimes burns, tans uniformly
    TYPE_IV = "type_iv"   # Burns minimally, always tans
    TYPE_V = "type_v"     # Very rarely burns, tans profusely
    TYPE_VI = "type_vi"   # Never burns, deeply pigmented

class MelanomaRisk(Enum):
    VERY_HIGH = "very_high"  # Multiple high-risk factors
    HIGH = "high"            # 2-3x population risk
    MODERATE = "moderate"    # 1.5-2x risk
    AVERAGE = "average"      # Population baseline
    LOW = "low"              # Protective factors

# =============================================================================
# MC1R - MASTER PIGMENTATION GENE
# =============================================================================

MC1R_MARKERS = {
    "rs1805007": {
        "gene": "MC1R",
        "variant": "R151C",
        "function": "Melanocortin 1 receptor - red hair/fair skin variant",
        "risk_allele": "T",
        "frequency": {"EUR": 0.10, "AFR": 0.001, "EAS": 0.001},
        "effect": {
            "TT": "R variant - Very high red hair/fair skin probability",
            "CT": "Carrier - Increased fair skin, freckling"
        },
        "category": "pigmentation",
        "phenotype": ["red_hair", "fair_skin", "freckling", "sunburn"],
        "evidence": "definitive",
        "pmid": ["8782809", "10854089", "18488028"],
        "melanoma_effect": "OR ~2.5 per allele",
        "actionable": {
            "TT": [
                "Very high probability of red/strawberry blonde hair",
                "Very fair skin - burns extremely easily",
                "HIGH melanoma risk (~6x)",
                "Sun protection essential - SPF 50+, protective clothing",
                "Regular skin cancer screening recommended",
                "Vitamin D synthesis may be adequate even with limited sun"
            ],
            "CT": [
                "Red hair carrier",
                "Increased freckling and fair skin tendency",
                "Elevated melanoma risk",
                "Good sun protection important"
            ]
        }
    },
    "rs1805008": {
        "gene": "MC1R",
        "variant": "R160W",
        "function": "MC1R high-penetrance variant",
        "risk_allele": "T",
        "frequency": {"EUR": 0.09},
        "effect": {
            "TT": "R variant - Red hair, very fair skin",
            "CT": "Carrier"
        },
        "category": "pigmentation",
        "phenotype": ["red_hair", "fair_skin"],
        "evidence": "definitive",
        "pmid": ["8782809", "18488028"],
        "melanoma_effect": "OR ~2.0 per allele"
    },
    "rs1805009": {
        "gene": "MC1R",
        "variant": "D294H",
        "function": "MC1R moderate effect variant",
        "risk_allele": "A",
        "frequency": {"EUR": 0.02},
        "effect": "Fair skin, melanoma risk",
        "category": "pigmentation",
        "evidence": "strong",
        "pmid": ["8782809"]
    },
    "rs1805005": {
        "gene": "MC1R",
        "variant": "V60L",
        "function": "MC1R low-penetrance (r) variant",
        "risk_allele": "T",
        "frequency": {"EUR": 0.15},
        "effect": "Mild effect on pigmentation",
        "category": "pigmentation",
        "evidence": "strong",
        "pmid": ["8782809", "18488028"],
        "melanoma_effect": "OR ~1.3 per allele"
    },
    "rs2228479": {
        "gene": "MC1R",
        "variant": "V92M",
        "function": "MC1R low-penetrance variant",
        "risk_allele": "A",
        "frequency": {"EUR": 0.10, "EAS": 0.30},
        "effect": "Mild effect on pigmentation",
        "category": "pigmentation",
        "evidence": "strong",
        "pmid": ["18488028"]
    },
    "rs11547464": {
        "gene": "MC1R",
        "variant": "R142H",
        "function": "MC1R variant",
        "risk_allele": "A",
        "frequency": {"EUR": 0.01},
        "effect": "Red hair variant",
        "category": "pigmentation",
        "evidence": "moderate",
        "pmid": ["8782809"]
    },
}

# =============================================================================
# OTHER PIGMENTATION GENES
# =============================================================================

PIGMENTATION_MARKERS = {
    **MC1R_MARKERS,
    "rs12913832": {
        "gene": "HERC2/OCA2",
        "variant": "Eye color major determinant",
        "function": "Regulates OCA2 expression - major eye color gene",
        "risk_allele": "G",  # G = brown, A = blue
        "frequency": {"EUR": 0.70, "AFR": 0.05, "EAS": 0.001},
        "effect": {
            "AA": "BLUE eyes ~75-85% probability",
            "AG": "Green/hazel likely, ~15% blue",
            "GG": "Brown eyes ~90%+ probability"
        },
        "category": "eye_color",
        "evidence": "definitive",
        "pmid": ["17952075", "18172690"],
        "actionable": {
            "AA": [
                "Blue eye color highly likely",
                "Blue eyes have less melanin - more UV sensitive",
                "Wear UV-protective sunglasses"
            ]
        }
    },
    "rs1800407": {
        "gene": "OCA2",
        "variant": "R419Q",
        "function": "P protein - melanin synthesis",
        "risk_allele": "T",  # A is blue modifier
        "frequency": {"EUR": 0.07},
        "effect": "Modifies toward blue/green eyes",
        "category": "eye_color",
        "evidence": "strong",
        "pmid": ["17236130", "18488028"]
    },
    "rs12203592": {
        "gene": "IRF4",
        "variant": "Interferon regulatory factor 4",
        "function": "Pigmentation regulation",
        "risk_allele": "T",
        "frequency": {"EUR": 0.16},
        "effect": {
            "TT": "Freckling, light eyes, sun sensitivity, early gray hair",
            "CT": "Moderate effect"
        },
        "category": "pigmentation",
        "evidence": "definitive",
        "pmid": ["18488028", "21378990"],
        "phenotypes": ["freckling", "light_eyes", "sun_sensitivity", "gray_hair"],
        "actionable": {
            "TT": [
                "High freckling tendency",
                "Light eye color more likely",
                "Earlier graying possible",
                "Good sun protection recommended"
            ]
        }
    },
    "rs1126809": {
        "gene": "TYR",
        "variant": "R402Q (Tyrosinase)",
        "function": "Rate-limiting enzyme in melanin synthesis",
        "risk_allele": "A",
        "frequency": {"EUR": 0.25},
        "effect": "Lighter pigmentation",
        "category": "pigmentation",
        "evidence": "strong",
        "pmid": ["12629594", "18488028"]
    },
    "rs16891982": {
        "gene": "SLC45A2",
        "variant": "L374F (MATP)",
        "function": "Membrane-associated transporter - melanin",
        "risk_allele": "C",  # European light skin allele
        "frequency": {"EUR": 0.95, "AFR": 0.05, "EAS": 0.05},
        "effect": {
            "CC": "European light skin variant",
            "GG": "Ancestral darker skin variant"
        },
        "category": "pigmentation",
        "evidence": "definitive",
        "pmid": ["15714523", "18488028"]
    },
    "rs1426654": {
        "gene": "SLC24A5",
        "variant": "A111T",
        "function": "Major European skin lightening variant",
        "risk_allele": "A",
        "frequency": {"EUR": 0.99, "AFR": 0.05, "EAS": 0.05, "SAS": 0.50},
        "effect": "Nearly fixed in Europeans - light skin",
        "category": "pigmentation",
        "evidence": "definitive",
        "pmid": ["16357253", "17182896"]
    },
}

# =============================================================================
# MELANOMA / SKIN CANCER
# =============================================================================

MELANOMA_MARKERS = {
    "rs910873": {
        "gene": "CDKN2A",
        "variant": "p16 variant",
        "function": "Cell cycle control - tumor suppressor",
        "risk_allele": "C",
        "frequency": {"EUR": 0.02},
        "effect": {
            "note": "Pathogenic mutations cause familial melanoma",
            "common_variant": "Associated with melanoma risk"
        },
        "category": "melanoma",
        "evidence": "strong",
        "pmid": ["8197451", "15928703"],
        "actionable": {
            "carrier": [
                "CDKN2A variants associated with familial melanoma",
                "Annual full-body skin exams recommended",
                "Self-exam monthly",
                "Strict sun protection"
            ]
        }
    },
    "rs401681": {
        "gene": "TERT-CLPTM1L",
        "variant": "5p15.33 region",
        "function": "Telomerase regulation",
        "risk_allele": "C",
        "frequency": {"EUR": 0.45},
        "effect": "Associated with melanoma and other cancers",
        "category": "melanoma",
        "evidence": "strong",
        "pmid": ["19578365", "18488028"]
    },
    # MC1R R variants above are major melanoma risk factors
}

# =============================================================================
# SKIN AGING
# =============================================================================

SKIN_AGING_MARKERS = {
    "rs1799750": {
        "gene": "MMP1",
        "variant": "1G/2G promoter",
        "function": "Matrix metalloproteinase 1 - collagen breakdown",
        "risk_allele": "G",  # 2G allele
        "frequency": {"EUR": 0.50},
        "effect": {
            "GG": "2G/2G - Higher MMP1 expression, faster collagen degradation",
            "note": "Accelerated skin aging, more wrinkles"
        },
        "category": "aging",
        "evidence": "moderate",
        "pmid": ["9843963", "11404010"],
        "actionable": {
            "GG": [
                "Higher collagen breakdown with UV exposure",
                "Sunscreen ESSENTIAL for aging prevention",
                "Consider retinoids for anti-aging",
                "Antioxidants (vitamin C, E) may help"
            ]
        }
    },
    "rs1800012": {
        "gene": "COL1A1",
        "variant": "Sp1 polymorphism",
        "function": "Collagen type I production",
        "risk_allele": "T",
        "frequency": {"EUR": 0.20},
        "effect": "Affects collagen production",
        "category": "aging",
        "evidence": "moderate",
        "pmid": ["16477526"]
    },
}

# =============================================================================
# ACNE, ROSACEA, PSORIASIS, ECZEMA
# =============================================================================

SKIN_CONDITION_MARKERS = {
    # Acne
    "rs2660753": {
        "gene": "1q25.3",
        "variant": "Acne GWAS hit",
        "function": "Unknown - acne susceptibility",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15},
        "effect": "Associated with severe acne",
        "category": "acne",
        "evidence": "moderate",
        "pmid": ["24927181"]
    },
    # Rosacea
    "rs763035": {
        "gene": "HLA region",
        "variant": "Rosacea susceptibility",
        "function": "Immune-related",
        "risk_allele": "A",
        "frequency": {"EUR": 0.20},
        "effect": "Associated with rosacea",
        "category": "rosacea",
        "evidence": "moderate",
        "pmid": ["25629664"]
    },
    # Psoriasis
    "rs10484554": {
        "gene": "HLA-C*06:02",
        "variant": "Major psoriasis risk",
        "function": "HLA class I - immune presentation",
        "risk_allele": "A",
        "frequency": {"EUR": 0.12},
        "effect": {
            "AA": "High psoriasis risk (HLA-Cw6)",
            "note": "Strongest genetic risk factor for psoriasis"
        },
        "category": "psoriasis",
        "evidence": "definitive",
        "pmid": ["17236132", "20953190"],
        "actionable": {
            "AA": [
                "Elevated psoriasis risk",
                "Early treatment if symptoms develop",
                "Biologics often effective (HLA-Cw6+ may respond better to IL-17 inhibitors)"
            ]
        }
    },
    # Eczema/Atopic dermatitis
    "rs61816761": {
        "gene": "FLG",
        "variant": "R501X (Filaggrin)",
        "function": "Skin barrier protein",
        "risk_allele": "A",
        "frequency": {"EUR": 0.04},
        "effect": {
            "AA": "Filaggrin null - very high eczema risk",
            "AG": "Carrier - elevated eczema risk"
        },
        "category": "eczema",
        "evidence": "definitive",
        "pmid": ["16550169", "19250312"],
        "actionable": {
            "carrier": [
                "Filaggrin deficiency - impaired skin barrier",
                "Higher eczema and dry skin risk",
                "Regular moisturizing essential",
                "Higher risk of ichthyosis vulgaris",
                "May also increase asthma risk"
            ]
        }
    },
    "rs2282458": {
        "gene": "FLG",
        "variant": "2282del4 proxy",
        "function": "Another filaggrin loss-of-function variant",
        "risk_allele": "deletion",
        "frequency": {"EUR": 0.02},
        "effect": "Filaggrin deficiency",
        "category": "eczema",
        "evidence": "definitive",
        "pmid": ["16550169"]
    },
}

# =============================================================================
# HAIR - COLOR, LOSS, GRAY, TEXTURE
# =============================================================================

HAIR_MARKERS = {
    # Hair color
    "rs12821256": {
        "gene": "KITLG",
        "variant": "Hair color modulator",
        "function": "Kit ligand - melanocyte development",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15},
        "effect": "Associated with blonde hair",
        "category": "hair_color",
        "evidence": "definitive",
        "pmid": ["18488028", "21378990"]
    },
    "rs28777": {
        "gene": "SLC45A2",
        "variant": "Hair color variant",
        "function": "Pigmentation transport",
        "risk_allele": "A",
        "frequency": {"EUR": 0.25},
        "effect": "Lighter hair color",
        "category": "hair_color",
        "evidence": "strong",
        "pmid": ["21378990"]
    },
    
    # Male pattern baldness
    "rs2180439": {
        "gene": "AR region (Xq12)",
        "variant": "Androgen receptor variant",
        "function": "Androgen sensitivity",
        "risk_allele": "T",
        "frequency": {"EUR": 0.50},
        "effect": {
            "T": "Associated with male pattern baldness",
            "note": "X-linked - males have one copy"
        },
        "category": "hair_loss",
        "evidence": "definitive",
        "pmid": ["18849991", "17963221"],
        "actionable": {
            "T_male": [
                "Elevated male pattern baldness risk",
                "Finasteride may be effective for prevention",
                "Minoxidil can help with regrowth",
                "Hair loss often from maternal grandfather"
            ]
        }
    },
    "rs6152": {
        "gene": "AR",
        "variant": "Androgen receptor (StuI)",
        "function": "Androgen receptor activity",
        "risk_allele": "G",
        "frequency": {"EUR": 0.40},
        "effect": "Male pattern baldness risk",
        "category": "hair_loss",
        "evidence": "strong",
        "pmid": ["17963221"]
    },
    "rs1160312": {
        "gene": "WNT10A",
        "variant": "Hair follicle development",
        "function": "WNT signaling in hair",
        "risk_allele": "A",
        "frequency": {"EUR": 0.40},
        "effect": "Associated with male pattern baldness",
        "category": "hair_loss",
        "evidence": "strong",
        "pmid": ["24174332"]
    },
    
    # Female hair loss
    "rs700518": {
        "gene": "CYP19A1",
        "variant": "Aromatase variant",
        "function": "Estrogen synthesis",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30},
        "effect": "Associated with female pattern hair loss",
        "category": "hair_loss",
        "evidence": "moderate",
        "pmid": ["18852991"]
    },
    
    # Gray hair timing
    "rs12203592": {
        "gene": "IRF4",
        "reference": "See PIGMENTATION_MARKERS",
        "gray_hair_effect": "TT associated with earlier graying"
    },
    "rs35264875": {
        "gene": "TERT",
        "variant": "Telomerase variant",
        "function": "Cellular aging",
        "risk_allele": "A",
        "frequency": {"EUR": 0.40},
        "effect": "Associated with premature graying",
        "category": "gray_hair",
        "evidence": "moderate",
        "pmid": ["21378990", "26926045"]
    },
    
    # Hair texture
    "rs3827760": {
        "gene": "EDAR",
        "variant": "370Val>Ala",
        "function": "Ectodysplasin A receptor - hair follicle",
        "risk_allele": "A",  # Derived Asian allele
        "frequency": {"EUR": 0.01, "EAS": 0.90, "AMR": 0.50},
        "effect": {
            "AA": "Thick, straight hair (East Asian type)",
            "note": "Also affects tooth shape and sweat glands"
        },
        "category": "hair_texture",
        "evidence": "definitive",
        "pmid": ["18561325", "24913232"],
        "actionable": {
            "AA": [
                "Thick, straight hair phenotype",
                "More active sweat glands",
                "Shovel-shaped incisors",
                "Common in East Asian and Native American ancestry"
            ]
        }
    },
    "rs11803731": {
        "gene": "TCHH",
        "variant": "Trichohyalin",
        "function": "Hair shaft structure",
        "risk_allele": "A",
        "frequency": {"EUR": 0.15},
        "effect": "Associated with straight vs curly hair",
        "category": "hair_texture",
        "evidence": "strong",
        "pmid": ["19578365", "19680544"]
    },
}

# =============================================================================
# COMBINE ALL SKIN/APPEARANCE MARKERS
# =============================================================================

SKIN_APPEARANCE_MARKERS = {
    **PIGMENTATION_MARKERS,
    **MELANOMA_MARKERS,
    **SKIN_AGING_MARKERS,
    **SKIN_CONDITION_MARKERS,
    **HAIR_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def calculate_skin_type(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Estimate Fitzpatrick skin type from genetics."""
    
    fair_score = 0
    
    # MC1R R variants (major effect)
    mc1r_r_variants = ["rs1805007", "rs1805008", "rs1805009"]
    for rs in mc1r_r_variants:
        geno = genotypes.get(rs, "CC")
        if geno[0] != "C":  # Has risk allele
            fair_score += 3 if geno == "TT" or geno == "AA" else 1.5
    
    # IRF4
    irf4 = genotypes.get("rs12203592", "CC")
    if irf4 == "TT":
        fair_score += 2
    elif irf4 == "CT":
        fair_score += 1
    
    # SLC45A2
    slc45a2 = genotypes.get("rs16891982", "CC")
    if slc45a2 == "CC":
        fair_score += 1  # European light variant
    elif slc45a2 == "GG":
        fair_score -= 2  # Ancestral dark variant
    
    # Determine skin type
    if fair_score >= 5:
        skin_type = SkinType.TYPE_I
        spf = "50+"
        description = "Very fair, always burns, never tans"
    elif fair_score >= 3:
        skin_type = SkinType.TYPE_II
        spf = "30-50"
        description = "Fair, usually burns, tans with difficulty"
    elif fair_score >= 1:
        skin_type = SkinType.TYPE_III
        spf = "15-30"
        description = "Medium, sometimes burns, tans uniformly"
    elif fair_score >= -1:
        skin_type = SkinType.TYPE_IV
        spf = "15"
        description = "Olive, rarely burns, always tans"
    elif fair_score >= -3:
        skin_type = SkinType.TYPE_V
        spf = "15"
        description = "Brown, very rarely burns, tans darkly"
    else:
        skin_type = SkinType.TYPE_VI
        spf = "15"
        description = "Dark brown/black, never burns"
    
    return {
        "skin_type": skin_type.value,
        "fitzpatrick_type": skin_type.name.replace("TYPE_", ""),
        "description": description,
        "recommended_spf": spf,
        "fair_score": round(fair_score, 1)
    }

def calculate_melanoma_risk(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate melanoma genetic risk."""
    
    risk_score = 0
    risk_factors = []
    
    # MC1R R variants (major risk)
    mc1r_count = 0
    for rs in ["rs1805007", "rs1805008", "rs1805009"]:
        geno = genotypes.get(rs, "CC")
        if geno[0] != "C":
            mc1r_count += 1 if geno[0] in "TA" else 0.5
    
    if mc1r_count >= 2:
        risk_score += 3
        risk_factors.append("Multiple MC1R R variants - high melanoma risk")
    elif mc1r_count >= 1:
        risk_score += 1.5
        risk_factors.append("MC1R R variant carrier")
    
    # IRF4
    irf4 = genotypes.get("rs12203592", "CC")
    if irf4 == "TT":
        risk_score += 1
        risk_factors.append("IRF4 TT - sun sensitivity and freckling")
    
    # CDKN2A
    cdkn2a = genotypes.get("rs910873", "TT")
    if cdkn2a in ["CC", "CT"]:
        risk_score += 1
        risk_factors.append("CDKN2A variant")
    
    # Determine risk level
    if risk_score >= 4:
        level = MelanomaRisk.VERY_HIGH
    elif risk_score >= 2.5:
        level = MelanomaRisk.HIGH
    elif risk_score >= 1.5:
        level = MelanomaRisk.MODERATE
    else:
        level = MelanomaRisk.AVERAGE
    
    return {
        "risk_level": level.value,
        "risk_factors": risk_factors,
        "recommendations": [
            "Annual full-body skin exam with dermatologist" if level in [MelanomaRisk.VERY_HIGH, MelanomaRisk.HIGH] else "Regular skin self-exams",
            "SPF 50+ daily, reapply every 2 hours in sun" if level != MelanomaRisk.AVERAGE else "SPF 30+ for sun exposure",
            "Avoid tanning beds - extremely high-risk" if level in [MelanomaRisk.VERY_HIGH, MelanomaRisk.HIGH] else None,
            "Protective clothing, seek shade 10am-4pm" if level != MelanomaRisk.AVERAGE else None,
            "Monthly self-exam using ABCDE criteria" if level != MelanomaRisk.AVERAGE else "Occasional self-exam"
        ],
        "pmid": ["18488028", "8782809"]
    }

def predict_eye_color(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Predict eye color from genetics."""
    
    herc2 = genotypes.get("rs12913832", "AG")
    oca2 = genotypes.get("rs1800407", "CC")
    
    if herc2 == "AA":
        primary = "blue"
        confidence = "high"
        probability = {"blue": 0.80, "green": 0.15, "brown": 0.05}
    elif herc2 == "AG":
        primary = "green/hazel"
        confidence = "moderate"
        probability = {"blue": 0.15, "green": 0.55, "brown": 0.30}
        if oca2 == "TT":
            primary = "green"
            probability = {"blue": 0.25, "green": 0.60, "brown": 0.15}
    else:  # GG
        primary = "brown"
        confidence = "high"
        probability = {"blue": 0.02, "green": 0.08, "brown": 0.90}
    
    return {
        "predicted_color": primary,
        "confidence": confidence,
        "probability_distribution": probability,
        "herc2_oca2_genotype": f"HERC2: {herc2}, OCA2: {oca2}",
        "pmid": ["17952075", "18172690"]
    }

def predict_hair_loss_risk(genotypes: Dict[str, str], sex: str = "male") -> Dict[str, Any]:
    """Predict hair loss risk based on genetics and sex."""
    
    if sex.lower() == "male":
        ar = genotypes.get("rs2180439", "C")
        ar2 = genotypes.get("rs6152", "A")
        wnt10a = genotypes.get("rs1160312", "G")
        
        risk_score = 0
        if ar == "T":
            risk_score += 2
        if ar2 == "G":
            risk_score += 1
        if wnt10a == "A":
            risk_score += 0.5
        
        if risk_score >= 3:
            risk = "HIGH"
            probability = "70-80%"
        elif risk_score >= 2:
            risk = "MODERATE"
            probability = "50-60%"
        else:
            risk = "LOW"
            probability = "<40%"
        
        return {
            "risk_level": risk,
            "probability_by_50": probability,
            "recommendations": [
                "Monitor hairline" if risk != "LOW" else None,
                "Finasteride may slow progression" if risk in ["HIGH", "MODERATE"] else None,
                "Minoxidil can help with regrowth" if risk in ["HIGH", "MODERATE"] else None
            ],
            "note": "Hair loss pattern often inherited from maternal grandfather"
        }
    else:  # female
        cyp19a1 = genotypes.get("rs700518", "GG")
        risk = "ELEVATED" if cyp19a1 == "AA" else "AVERAGE"
        
        return {
            "risk_level": risk,
            "recommendations": [
                "Check thyroid and iron levels if experiencing loss",
                "Minoxidil (2%) approved for female pattern loss"
            ]
        }

def generate_skin_appearance_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate comprehensive skin and appearance report."""
    
    skin_type = calculate_skin_type(genotypes)
    melanoma_risk = calculate_melanoma_risk(genotypes)
    eye_color = predict_eye_color(genotypes)
    
    return {
        "skin": skin_type,
        "melanoma_risk": melanoma_risk,
        "eye_color": eye_color,
        "traits": {
            "likely_freckles": genotypes.get("rs12203592", "CC") == "TT",
            "red_hair_carrier": any(genotypes.get(rs, "CC")[0] != "C" 
                                   for rs in ["rs1805007", "rs1805008"]),
        },
        "markers_analyzed": sum(1 for rs in SKIN_APPEARANCE_MARKERS if rs in genotypes),
    }

# Export
__all__ = [
    'SKIN_APPEARANCE_MARKERS',
    'PIGMENTATION_MARKERS',
    'MC1R_MARKERS',
    'MELANOMA_MARKERS',
    'SKIN_AGING_MARKERS',
    'SKIN_CONDITION_MARKERS',
    'HAIR_MARKERS',
    'SkinType',
    'MelanomaRisk',
    'calculate_skin_type',
    'calculate_melanoma_risk',
    'predict_eye_color',
    'predict_hair_loss_risk',
    'generate_skin_appearance_report',
]
