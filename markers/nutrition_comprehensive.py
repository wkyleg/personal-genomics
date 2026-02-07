"""
Comprehensive Nutrition Genetics v5.0

Complete coverage of:
- Fat metabolism (APOE, APOA2, PPARG, FABP2)
- Carb tolerance (TCF7L2, IRS1, PPARG)
- Sugar/insulin response
- Sodium sensitivity
- Histamine intolerance (DAO, HNMT, MAO)
- Oxalate metabolism
- Gluten/Celiac (complete HLA-DQ panel)
- Lactose intolerance
- Fructose malabsorption
- Alcohol metabolism
- Taste genetics (bitter, sweet, fat, umami, salt)
- Caffeine detailed

All markers with PMID references and dietary recommendations.
"""

from enum import Enum
from typing import Dict, List, Any, Optional

class ToleranceLevel(Enum):
    EXCELLENT = "excellent"
    GOOD = "good"
    MODERATE = "moderate"
    POOR = "poor"
    VERY_POOR = "very_poor"

class DietaryRecommendation(Enum):
    LOW_FAT = "low_fat"
    MEDITERRANEAN = "mediterranean"
    LOW_CARB = "low_carb"
    BALANCED = "balanced"
    KETO_COMPATIBLE = "keto_compatible"

# =============================================================================
# FAT METABOLISM - APOE
# =============================================================================

APOE_MARKERS = {
    "rs429358": {
        "gene": "APOE",
        "variant": "APOE ε4 determinant (C112R)",
        "function": "Determines ε4 vs ε2/ε3",
        "risk_allele": "C",  # C = ε4
        "frequency": {"EUR": 0.15, "AFR": 0.25, "EAS": 0.08},
        "effect": "Combined with rs7412 determines APOE genotype",
        "category": "fat_metabolism",
        "evidence": "definitive",
        "pmid": ["8346443", "19734902", "25849268"],
        "note": "Must combine with rs7412 for full APOE typing"
    },
    "rs7412": {
        "gene": "APOE",
        "variant": "APOE ε2 determinant (R158C)",
        "function": "Determines ε2 vs ε3/ε4",
        "risk_allele": "T",  # T = ε2
        "frequency": {"EUR": 0.08, "AFR": 0.11, "EAS": 0.05},
        "effect": "Combined with rs429358 determines APOE genotype",
        "category": "fat_metabolism",
        "evidence": "definitive",
        "pmid": ["8346443"],
        "apoe_determination": {
            # rs429358/rs7412 combinations
            "CC/TT": "ε2/ε4 (rare, conflicting)",
            "CC/TC": "ε3/ε4",
            "CC/CC": "ε4/ε4",
            "CT/TT": "ε2/ε2",
            "CT/TC": "ε2/ε3",
            "CT/CC": "ε2/ε4",
            "TT/TT": "ε2/ε2",
            "TT/TC": "ε2/ε3",
            "TT/CC": "ε3/ε3"
        }
    },
}

APOE_DIET_RECOMMENDATIONS = {
    "e2/e2": {
        "saturated_fat": ToleranceLevel.GOOD,
        "carb_tolerance": ToleranceLevel.POOR,  # Hyperlipoproteinemia type III risk
        "diet_type": DietaryRecommendation.BALANCED,
        "recommendations": [
            "May tolerate saturated fat better than other genotypes",
            "Watch for elevated triglycerides with high carb",
            "Moderate carbohydrate intake recommended",
            "Generally favorable cardiovascular risk"
        ],
        "pmid": ["25849268"]
    },
    "e2/e3": {
        "saturated_fat": ToleranceLevel.GOOD,
        "carb_tolerance": ToleranceLevel.MODERATE,
        "diet_type": DietaryRecommendation.BALANCED,
        "recommendations": [
            "Flexible dietary response",
            "Can tolerate moderate saturated fat",
            "Lower LDL than average"
        ]
    },
    "e3/e3": {
        "saturated_fat": ToleranceLevel.MODERATE,
        "carb_tolerance": ToleranceLevel.MODERATE,
        "diet_type": DietaryRecommendation.BALANCED,
        "recommendations": [
            "Most common genotype - average response",
            "Standard dietary recommendations apply",
            "Mediterranean diet well-suited"
        ]
    },
    "e3/e4": {
        "saturated_fat": ToleranceLevel.POOR,
        "carb_tolerance": ToleranceLevel.GOOD,
        "diet_type": DietaryRecommendation.LOW_FAT,
        "recommendations": [
            "LIMIT saturated fat - more responsive to dietary fat",
            "Lower saturated fat (<7% calories) recommended",
            "Mediterranean diet preferable",
            "May tolerate carbohydrates better than fat",
            "Higher LDL response to saturated fat"
        ],
        "pmid": ["10479231", "21139125"]
    },
    "e4/e4": {
        "saturated_fat": ToleranceLevel.VERY_POOR,
        "carb_tolerance": ToleranceLevel.GOOD,
        "diet_type": DietaryRecommendation.LOW_FAT,
        "recommendations": [
            "STRICT saturated fat restriction essential",
            "Highest LDL response to saturated fat",
            "Mediterranean diet strongly recommended",
            "Consider plant-based emphasis",
            "Omega-3 fatty acids especially important",
            "Alzheimer's risk - consider brain-healthy diet (MIND diet)"
        ],
        "pmid": ["10479231", "21139125", "25849268"]
    },
    "e2/e4": {
        "saturated_fat": ToleranceLevel.MODERATE,
        "carb_tolerance": ToleranceLevel.MODERATE,
        "diet_type": DietaryRecommendation.MEDITERRANEAN,
        "recommendations": [
            "Mixed signals - moderate all macros",
            "Mediterranean approach is safest"
        ]
    }
}

FAT_METABOLISM_MARKERS = {
    **APOE_MARKERS,
    "rs5082": {
        "gene": "APOA2",
        "variant": "-265T>C promoter",
        "function": "ApoA-II expression - saturated fat response",
        "risk_allele": "C",
        "frequency": {"EUR": 0.25, "EAS": 0.10, "AFR": 0.35},
        "effect": {
            "CC": "Saturated fat causes MORE weight gain",
            "CT": "Moderate effect",
            "TT": "Resistant to saturated fat-induced weight gain"
        },
        "category": "fat_metabolism",
        "evidence": "strong",
        "pmid": ["19662675", "19424078"],
        "actionable": {
            "CC": [
                "HIGH sensitivity to saturated fat → weight gain",
                "Limit saturated fat to <22g/day (<10% calories)",
                "Replace with unsaturated fats (olive oil, nuts, fish)",
                "Coconut oil, butter, red meat → weight gain for you"
            ],
            "TT": [
                "Less affected by saturated fat for weight",
                "Still consider cardiovascular effects"
            ]
        }
    },
    "rs1801282": {
        "gene": "PPARG",
        "variant": "Pro12Ala",
        "function": "Peroxisome proliferator-activated receptor gamma",
        "risk_allele": "C",  # Pro allele = C
        "frequency": {"EUR": 0.88, "EAS": 0.95, "AFR": 0.99},
        "effect": {
            "CC": "Pro/Pro - Higher T2D risk, insulin resistance",
            "CG": "Pro/Ala - Intermediate",
            "GG": "Ala/Ala - Protective against T2D, better insulin sensitivity"
        },
        "category": "fat_metabolism",
        "evidence": "definitive",
        "pmid": ["9614613", "11574435", "19584936"],
        "actionable": {
            "CC": [
                "Higher T2D and insulin resistance risk",
                "Weight management especially important",
                "May benefit from lower fat diet",
                "Responds well to monounsaturated fat (MUFA)"
            ],
            "GG": [
                "Ala/Ala protective against T2D",
                "Better insulin sensitivity",
                "Can tolerate higher fat diets"
            ]
        }
    },
    "rs1799883": {
        "gene": "FABP2",
        "variant": "Ala54Thr",
        "function": "Fatty acid binding protein 2 - intestinal fat absorption",
        "risk_allele": "A",  # Thr allele
        "frequency": {"EUR": 0.28, "EAS": 0.35, "AFR": 0.20},
        "effect": {
            "AA": "Thr/Thr - Higher fat absorption, higher insulin response",
            "AG": "Intermediate",
            "GG": "Ala/Ala - Normal fat absorption"
        },
        "category": "fat_metabolism",
        "evidence": "strong",
        "pmid": ["7550340", "18772348"],
        "actionable": {
            "AA": [
                "Higher intestinal fat absorption",
                "May benefit from lower dietary fat",
                "Higher postprandial triglycerides"
            ]
        }
    },
}

# =============================================================================
# CARBOHYDRATE / INSULIN RESPONSE
# =============================================================================

CARB_MARKERS = {
    "rs7903146": {
        "gene": "TCF7L2",
        "variant": "TCF7L2 T2D risk variant",
        "function": "Transcription factor affecting insulin secretion",
        "risk_allele": "T",
        "frequency": {"EUR": 0.30, "AFR": 0.30, "EAS": 0.05},
        "effect": {
            "TT": "1.8x T2D risk - impaired insulin secretion",
            "CT": "1.4x T2D risk",
            "CC": "Normal risk"
        },
        "category": "carb_tolerance",
        "evidence": "definitive",
        "pmid": ["16415884", "17463248", "19174826"],
        "actionable": {
            "TT": [
                "HIGHEST T2D genetic risk (strongest known variant)",
                "Limit refined carbohydrates and added sugars",
                "Focus on low glycemic index foods",
                "Regular blood glucose monitoring recommended",
                "Physical activity especially important"
            ],
            "CT": [
                "Elevated T2D risk",
                "Moderate carbohydrate restriction beneficial"
            ]
        }
    },
    "rs2943641": {
        "gene": "IRS1",
        "variant": "Insulin receptor substrate 1",
        "function": "Insulin signaling pathway",
        "risk_allele": "C",
        "frequency": {"EUR": 0.60, "EAS": 0.70, "AFR": 0.45},
        "effect": {
            "CC": "Insulin resistance tendency",
            "TT": "Better insulin sensitivity"
        },
        "category": "carb_tolerance",
        "evidence": "strong",
        "pmid": ["22231480", "18556340"],
        "actionable": {
            "CC": [
                "Tendency toward insulin resistance",
                "Lower glycemic load diet beneficial",
                "Exercise improves insulin sensitivity"
            ]
        }
    },
    "rs13266634": {
        "gene": "SLC30A8",
        "variant": "R325W (zinc transporter)",
        "function": "Zinc transport in pancreatic beta cells",
        "risk_allele": "C",  # R allele = C
        "frequency": {"EUR": 0.70, "EAS": 0.60, "AFR": 0.95},
        "effect": {
            "CC": "R/R - Higher T2D risk",
            "TT": "W/W - Protective against T2D (paradoxically loss of function)"
        },
        "category": "carb_tolerance",
        "evidence": "definitive",
        "pmid": ["17463249", "24136000"]
    },
    "rs5219": {
        "gene": "KCNJ11",
        "variant": "E23K",
        "function": "Potassium channel in pancreatic beta cells",
        "risk_allele": "T",  # K allele
        "frequency": {"EUR": 0.35, "EAS": 0.40, "AFR": 0.10},
        "effect": {
            "TT": "K/K - Reduced insulin secretion, T2D risk",
            "CC": "E/E - Normal insulin secretion"
        },
        "category": "carb_tolerance",
        "evidence": "strong",
        "pmid": ["15677319", "17463248"]
    },
}

# =============================================================================
# SODIUM SENSITIVITY
# =============================================================================

SODIUM_MARKERS = {
    "rs4340": {
        "gene": "ACE",
        "variant": "I/D polymorphism",
        "function": "Angiotensin-converting enzyme",
        "risk_allele": "D",
        "frequency": {"EUR": 0.52, "AFR": 0.60},
        "effect": {
            "DD": "Salt-sensitive blood pressure",
            "II": "Less salt-sensitive"
        },
        "category": "sodium",
        "evidence": "strong",
        "pmid": ["8675673", "16567525"],
        "actionable": {
            "DD": [
                "Salt-sensitive blood pressure",
                "Limit sodium to <2000mg/day",
                "DASH diet particularly effective"
            ]
        }
    },
    "rs4961": {
        "gene": "ADD1",
        "variant": "Gly460Trp (alpha-adducin)",
        "function": "Sodium reabsorption in kidney",
        "risk_allele": "T",  # Trp allele
        "frequency": {"EUR": 0.22, "EAS": 0.50},
        "effect": {
            "TT": "Trp/Trp - HIGHLY salt-sensitive BP",
            "GT": "Moderate salt sensitivity",
            "GG": "Normal salt response"
        },
        "category": "sodium",
        "evidence": "strong",
        "pmid": ["10591398", "16567525"],
        "actionable": {
            "TT": [
                "HIGHLY salt-sensitive",
                "STRICT sodium restriction (<1500mg/day)",
                "Diuretics particularly effective",
                "Potassium-rich foods beneficial"
            ]
        }
    },
    "rs699": {
        "gene": "AGT",
        "variant": "M235T",
        "function": "Angiotensinogen levels",
        "risk_allele": "G",  # T allele
        "frequency": {"EUR": 0.42, "AFR": 0.90},
        "effect": {
            "GG": "Higher AGT, salt-sensitive",
            "AA": "Lower AGT"
        },
        "category": "sodium",
        "evidence": "strong",
        "pmid": ["1641269"]
    },
}

# =============================================================================
# HISTAMINE INTOLERANCE
# =============================================================================

HISTAMINE_MARKERS = {
    "rs1049793": {
        "gene": "AOC1/DAO",
        "variant": "Diamine oxidase C/G",
        "function": "Primary enzyme degrading histamine in gut",
        "risk_allele": "G",
        "frequency": {"EUR": 0.25, "EAS": 0.15},
        "effect": {
            "GG": "Reduced DAO activity - histamine intolerance likely",
            "CG": "Intermediate DAO activity"
        },
        "category": "histamine",
        "evidence": "moderate",
        "pmid": ["23579755", "28158556"],
        "actionable": {
            "GG": [
                "HIGH histamine intolerance risk",
                "Avoid high-histamine foods (aged cheese, wine, fermented foods)",
                "Fresh foods preferred over leftovers",
                "DAO enzyme supplements may help before meals",
                "Symptoms: headaches, flushing, GI upset, hives"
            ]
        }
    },
    "rs1049742": {
        "gene": "AOC1/DAO",
        "variant": "DAO His645Asp",
        "function": "DAO enzyme activity",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15},
        "effect": {
            "CC": "Reduced DAO activity"
        },
        "category": "histamine",
        "evidence": "moderate",
        "pmid": ["23579755"]
    },
    "rs11558538": {
        "gene": "HNMT",
        "variant": "Thr105Ile",
        "function": "Histamine N-methyltransferase - intracellular histamine degradation",
        "risk_allele": "A",  # Ile allele
        "frequency": {"EUR": 0.10, "EAS": 0.15},
        "effect": {
            "AA": "Ile/Ile - Reduced HNMT activity, histamine buildup"
        },
        "category": "histamine",
        "evidence": "moderate",
        "pmid": ["9734076", "28158556"],
        "actionable": {
            "AA": [
                "Reduced intracellular histamine breakdown",
                "Combined with DAO variants increases intolerance risk",
                "May affect response to antihistamines"
            ]
        }
    },
    "rs1137070": {
        "gene": "MAOB",
        "variant": "MAO-B intron 13",
        "function": "Monoamine oxidase B - tyramine/histamine",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35},
        "effect": "Affects tyramine and histamine metabolism",
        "category": "histamine",
        "evidence": "moderate",
        "pmid": ["8105684"]
    },
}

# =============================================================================
# OXALATE METABOLISM
# =============================================================================

OXALATE_MARKERS = {
    "rs34116584": {
        "gene": "AGXT",
        "variant": "Alanine-glyoxylate aminotransferase",
        "function": "Oxalate synthesis pathway - primary hyperoxaluria gene",
        "risk_allele": "A",
        "frequency": {"EUR": 0.10},
        "effect": {
            "note": "Pathogenic mutations cause primary hyperoxaluria",
            "common_variants": "Common variants may affect oxalate levels"
        },
        "category": "oxalate",
        "evidence": "moderate",
        "pmid": ["10958762"],
        "actionable": {
            "risk": [
                "If kidney stone history, consider low-oxalate diet",
                "High oxalate: spinach, rhubarb, beets, nuts, chocolate",
                "Adequate calcium intake helps bind oxalate",
                "Hydration critical"
            ]
        }
    },
    "rs11677394": {
        "gene": "GRHPR",
        "variant": "Glyoxylate reductase",
        "function": "Alternative oxalate pathway",
        "risk_allele": "T",
        "frequency": {"EUR": 0.15},
        "effect": "Associated with oxalate metabolism",
        "category": "oxalate",
        "evidence": "moderate",
        "pmid": ["15210650"]
    },
}

# =============================================================================
# GLUTEN / CELIAC - COMPLETE HLA PANEL
# =============================================================================

CELIAC_MARKERS = {
    "rs2187668": {
        "gene": "HLA-DQA1*05",
        "variant": "DQ2.5 alpha chain",
        "function": "HLA-DQ2.5 component (95% of celiac)",
        "risk_allele": "T",
        "frequency": {"EUR": 0.20, "AFR": 0.10, "EAS": 0.05},
        "effect": {
            "TT": "Homozygous DQ2.5 component - HIGHEST celiac risk",
            "CT": "Heterozygous - elevated risk"
        },
        "category": "gluten",
        "evidence": "definitive",
        "pmid": ["18509540", "20190752", "24867074"],
        "actionable": {
            "TT": [
                "HIGH celiac disease risk (DQ2.5 likely)",
                "Screen with tTG-IgA if any GI symptoms",
                "Do NOT start gluten-free diet before testing",
                "First-degree relatives should be screened"
            ]
        }
    },
    "rs7454108": {
        "gene": "HLA-DQB1*02",
        "variant": "DQ2.5 beta chain",
        "function": "Second component of DQ2.5",
        "risk_allele": "C",
        "frequency": {"EUR": 0.25},
        "effect": "Combined with DQA1*05 for full DQ2.5",
        "category": "gluten",
        "evidence": "definitive",
        "pmid": ["18509540"]
    },
    "rs7775228": {
        "gene": "HLA-DQB1*03:02",
        "variant": "DQ8 marker",
        "function": "HLA-DQ8 (5-10% of celiac without DQ2)",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15, "AMR": 0.25},
        "effect": {
            "CC": "DQ8 positive - celiac risk if no DQ2",
            "note": "DQ8 sufficient for celiac without DQ2"
        },
        "category": "gluten",
        "evidence": "definitive",
        "pmid": ["18509540", "20190752"],
        "actionable": {
            "CC": [
                "DQ8 positive - celiac possible even without DQ2",
                "5-10% of celiac patients are DQ8+/DQ2-",
                "Screen if symptomatic"
            ]
        }
    },
    "rs2395182": {
        "gene": "HLA-DQA1*02",
        "variant": "DQ2.2 alpha component",
        "function": "Part of DQ2.2 haplotype",
        "risk_allele": "T",
        "frequency": {"EUR": 0.25},
        "effect": "DQ2.2 lower risk than DQ2.5 but still elevated",
        "category": "gluten",
        "evidence": "strong",
        "pmid": ["24867074"]
    },
}

# =============================================================================
# LACTOSE INTOLERANCE
# =============================================================================

LACTOSE_MARKERS = {
    "rs4988235": {
        "gene": "MCM6/LCT",
        "variant": "-13910 C>T (European persistence)",
        "function": "Lactase persistence regulatory variant",
        "risk_allele": "C",
        "frequency": {"EUR": 0.25, "AFR": 0.85, "EAS": 0.95},
        "effect": {
            "CC": "Lactase non-persistent - lactose intolerant as adult",
            "CT": "Heterozygous - usually lactase persistent",
            "TT": "Lactase persistent - can digest lactose lifelong"
        },
        "category": "lactose",
        "evidence": "definitive",
        "pmid": ["11788828", "12019640"],
        "actionable": {
            "CC": [
                "Lactose intolerant (primary hypolactasia)",
                "Avoid or limit dairy, or use lactase supplements",
                "Hard cheeses and yogurt usually tolerated (low lactose)",
                "Ensure adequate calcium from non-dairy sources"
            ],
            "TT": [
                "Lactase persistent - can digest milk throughout life",
                "Dairy is appropriate calcium source"
            ]
        }
    },
    "rs182549": {
        "gene": "MCM6/LCT",
        "variant": "-22018 G>A",
        "function": "Secondary lactase persistence marker",
        "risk_allele": "C",
        "frequency": {"EUR": 0.25},
        "effect": "In LD with rs4988235",
        "category": "lactose",
        "evidence": "strong",
        "pmid": ["12019640"]
    },
    "rs145946881": {
        "gene": "LCT",
        "variant": "-14010 G>C (East African persistence)",
        "function": "Independent lactase persistence in East Africa",
        "risk_allele": "C",
        "frequency": {"AFR_EA": 0.15},
        "effect": "Lactase persistence in some East African populations",
        "category": "lactose",
        "evidence": "definitive",
        "pmid": ["17159977"]
    },
}

# =============================================================================
# FRUCTOSE MALABSORPTION
# =============================================================================

FRUCTOSE_MARKERS = {
    "rs76917243": {
        "gene": "ALDOB",
        "variant": "Aldolase B (A149P)",
        "function": "Fructose metabolism - HFI gene",
        "risk_allele": "C",
        "frequency": {"EUR": 0.005},  # Rare but important
        "effect": {
            "note": "Hereditary fructose intolerance (HFI) mutations",
            "CC": "If pathogenic variant - avoid fructose entirely"
        },
        "category": "fructose",
        "evidence": "definitive",
        "pmid": ["21907358"],
        "actionable": {
            "pathogenic": [
                "Hereditary fructose intolerance - STRICT fructose avoidance",
                "Avoid fructose, sucrose, sorbitol",
                "Can cause hypoglycemia, liver/kidney damage"
            ]
        }
    },
    "rs1501299": {
        "gene": "ADIPOQ",
        "variant": "Adiponectin (affects fructose response)",
        "function": "Adiponectin levels - metabolic regulation",
        "risk_allele": "G",
        "frequency": {"EUR": 0.30, "EAS": 0.45},
        "effect": "May affect fructose-induced metabolic effects",
        "category": "fructose",
        "evidence": "moderate",
        "pmid": ["15731354"]
    },
}

# =============================================================================
# ALCOHOL METABOLISM
# =============================================================================

ALCOHOL_MARKERS = {
    "rs1229984": {
        "gene": "ADH1B",
        "variant": "Arg48His (*2)",
        "function": "Alcohol dehydrogenase - ethanol to acetaldehyde",
        "risk_allele": "T",  # His allele = fast metabolism
        "frequency": {"EUR": 0.05, "EAS": 0.70, "AFR": 0.02},
        "effect": {
            "TT": "His/His - VERY fast alcohol metabolism → acetaldehyde buildup",
            "CT": "Fast metabolism - flushing, nausea from alcohol",
            "CC": "Normal/slow alcohol metabolism"
        },
        "category": "alcohol",
        "evidence": "definitive",
        "pmid": ["7635462", "19124506"],
        "actionable": {
            "TT": [
                "Extreme Asian flush response likely",
                "Rapid acetaldehyde buildup - unpleasant reaction",
                "PROTECTIVE against alcoholism",
                "Still at risk for acetaldehyde-related cancer with drinking"
            ],
            "CT": [
                "Flushing reaction to alcohol",
                "Lower alcoholism risk"
            ]
        }
    },
    "rs698": {
        "gene": "ADH1C",
        "variant": "Ile349Val (*2)",
        "function": "ADH1C enzyme speed",
        "risk_allele": "T",  # Val = slower
        "frequency": {"EUR": 0.40, "EAS": 0.10},
        "effect": {
            "AA": "Ile/Ile - Fast ADH1C",
            "TT": "Val/Val - Slow ADH1C"
        },
        "category": "alcohol",
        "evidence": "moderate",
        "pmid": ["7635462"]
    },
    "rs671": {
        "gene": "ALDH2",
        "variant": "*2 (Glu504Lys)",
        "function": "Aldehyde dehydrogenase - clears acetaldehyde",
        "risk_allele": "A",  # Lys = inactive
        "frequency": {"EUR": 0.001, "EAS": 0.25, "AFR": 0.001},
        "effect": {
            "AA": "*2/*2 - NO ALDH2 activity - SEVERE flush, can't drink",
            "GA": "*1/*2 - Reduced activity - moderate flushing, cancer risk",
            "GG": "*1/*1 - Normal ALDH2"
        },
        "category": "alcohol",
        "evidence": "definitive",
        "pmid": ["1180215", "19124506", "26876464"],
        "actionable": {
            "AA": [
                "CANNOT metabolize acetaldehyde",
                "Severe flushing, nausea, racing heart with any alcohol",
                "AVOID alcohol entirely",
                "Using alcohol anyway increases esophageal cancer risk 10x"
            ],
            "GA": [
                "Reduced ALDH2 - flushing with alcohol",
                "If you drink despite flushing: HIGH esophageal cancer risk",
                "Consider avoiding alcohol"
            ]
        }
    },
}

# =============================================================================
# TASTE GENETICS
# =============================================================================

TASTE_MARKERS = {
    # Bitter taste
    "rs713598": {
        "gene": "TAS2R38",
        "variant": "A49P (PAV/AVI)",
        "function": "Bitter taste receptor - PTC/PROP tasting",
        "risk_allele": "G",  # AVI (non-taster) = G
        "frequency": {"EUR": 0.45, "AFR": 0.15, "EAS": 0.30},
        "effect": {
            "CC": "PAV/PAV - Super-taster (very sensitive to bitter)",
            "CG": "PAV/AVI - Medium taster",
            "GG": "AVI/AVI - Non-taster (can't taste PTC/PROP)"
        },
        "category": "taste",
        "evidence": "definitive",
        "pmid": ["12595690", "17570398"],
        "actionable": {
            "CC": [
                "Super-taster - very sensitive to bitter compounds",
                "May dislike: Brussels sprouts, kale, grapefruit, coffee, beer",
                "May need to mask bitter taste (fat, salt, sweetness) in healthy foods"
            ],
            "GG": [
                "Non-taster - can't detect PTC/PROP bitterness",
                "May enjoy bitter vegetables more easily",
                "May be less sensitive to caffeine bitterness"
            ]
        }
    },
    "rs1726866": {
        "gene": "TAS2R38",
        "variant": "V262A",
        "function": "Second TAS2R38 variant",
        "risk_allele": "G",
        "frequency": {"EUR": 0.45},
        "effect": "Part of PAV/AVI haplotype",
        "category": "taste",
        "evidence": "definitive",
        "pmid": ["12595690"]
    },
    "rs10246939": {
        "gene": "TAS2R38",
        "variant": "I296V",
        "function": "Third TAS2R38 variant",
        "risk_allele": "G",
        "frequency": {"EUR": 0.45},
        "effect": "Completes PAV/AVI determination",
        "category": "taste",
        "evidence": "definitive",
        "pmid": ["12595690"]
    },
    
    # Sweet taste
    "rs35874116": {
        "gene": "TAS1R2",
        "variant": "Sweet receptor variant",
        "function": "Sweet taste receptor",
        "risk_allele": "T",
        "frequency": {"EUR": 0.20, "AFR": 0.40},
        "effect": {
            "TT": "Reduced sweet sensitivity - may crave more sugar"
        },
        "category": "taste",
        "evidence": "moderate",
        "pmid": ["17095015", "21487388"],
        "actionable": {
            "TT": [
                "Lower sweet sensitivity",
                "May add more sugar to foods",
                "Be mindful of sugar intake"
            ]
        }
    },
    "rs307355": {
        "gene": "TAS1R3",
        "variant": "Sweet/Umami receptor",
        "function": "Sweet and umami detection",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35},
        "effect": "Affects sweet and umami perception",
        "category": "taste",
        "evidence": "moderate",
        "pmid": ["17095015"]
    },
    
    # Fat taste
    "rs1761667": {
        "gene": "CD36",
        "variant": "Fat taste receptor",
        "function": "Fatty acid detection on tongue",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35, "AFR": 0.60},
        "effect": {
            "AA": "Reduced fat taste sensitivity - may consume more fat",
            "GG": "High fat taste sensitivity"
        },
        "category": "taste",
        "evidence": "strong",
        "pmid": ["22570252", "26232584"],
        "actionable": {
            "AA": [
                "Lower fat taste sensitivity",
                "May prefer higher-fat foods",
                "Be mindful of fat intake - may overconsume",
                "Associated with higher BMI in some studies"
            ],
            "GG": [
                "Sensitive to fat taste",
                "May naturally prefer lower-fat foods"
            ]
        }
    },
    
    # Umami
    "rs34160967": {
        "gene": "TAS1R1",
        "variant": "Umami receptor",
        "function": "Glutamate (umami) taste",
        "risk_allele": "A",
        "frequency": {"EUR": 0.15, "EAS": 0.40},
        "effect": "Affects umami sensitivity",
        "category": "taste",
        "evidence": "moderate",
        "pmid": ["16051168"]
    },
    
    # Salt
    "rs239345": {
        "gene": "TRPV1",
        "variant": "Capsaicin/salt receptor variant",
        "function": "Affects salt taste sensitivity",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30},
        "effect": "Associated with salt taste perception",
        "category": "taste",
        "evidence": "limited",
        "pmid": ["23497217"]
    },
    
    # Cilantro
    "rs72921001": {
        "gene": "OR6A2",
        "variant": "Olfactory receptor",
        "function": "Detects aldehydes in cilantro",
        "risk_allele": "A",
        "frequency": {"EUR": 0.15, "EAS": 0.10, "SAS": 0.05},
        "effect": {
            "AA": "Cilantro tastes like soap - strong aversion",
            "GA": "May have mild soapy taste from cilantro",
            "GG": "Enjoys cilantro normally"
        },
        "category": "taste",
        "evidence": "moderate",
        "pmid": ["22927850"],
        "actionable": {
            "AA": [
                "Cilantro tastes like soap (genetic, not preference)",
                "Use parsley or other herbs as substitute",
                "Not a 'picky eater' - it's genetic!"
            ]
        }
    },
}

# =============================================================================
# COMBINE ALL NUTRITION MARKERS
# =============================================================================

NUTRITION_COMPREHENSIVE_MARKERS = {
    **FAT_METABOLISM_MARKERS,
    **CARB_MARKERS,
    **SODIUM_MARKERS,
    **HISTAMINE_MARKERS,
    **OXALATE_MARKERS,
    **CELIAC_MARKERS,
    **LACTOSE_MARKERS,
    **FRUCTOSE_MARKERS,
    **ALCOHOL_MARKERS,
    **TASTE_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def determine_apoe_genotype(genotypes: Dict[str, str]) -> str:
    """Determine APOE genotype from rs429358 and rs7412."""
    rs429358 = genotypes.get("rs429358", "TT")
    rs7412 = genotypes.get("rs7412", "CC")
    
    # Mapping
    combinations = {
        ("TT", "TT"): "ε2/ε2",
        ("TT", "CT"): "ε2/ε3",
        ("TT", "CC"): "ε3/ε3",
        ("CT", "TT"): "ε2/ε2",
        ("CT", "CT"): "ε2/ε3",
        ("CT", "CC"): "ε2/ε4",
        ("CC", "TT"): "ε2/ε4",
        ("CC", "CT"): "ε3/ε4",
        ("CC", "CC"): "ε4/ε4",
    }
    
    return combinations.get((rs429358, rs7412), "ε3/ε3")

def generate_nutrition_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate comprehensive nutrition genetics report."""
    
    # APOE
    apoe = determine_apoe_genotype(genotypes)
    apoe_recs = APOE_DIET_RECOMMENDATIONS.get(apoe.replace("ε", "e"), APOE_DIET_RECOMMENDATIONS["e3/e3"])
    
    # Saturated fat
    apoa2 = genotypes.get("rs5082", "TC")
    sat_fat_sensitive = apoa2 == "CC" or "e4" in apoe
    
    # Carb tolerance
    tcf7l2 = genotypes.get("rs7903146", "CC")
    carb_tolerance = ToleranceLevel.POOR if tcf7l2 == "TT" else ToleranceLevel.MODERATE if tcf7l2 == "CT" else ToleranceLevel.GOOD
    
    # Salt sensitivity
    add1 = genotypes.get("rs4961", "GG")
    ace = genotypes.get("rs4340", "ID")
    salt_sensitive = add1 == "TT" or ace == "DD"
    
    # Lactose
    lct = genotypes.get("rs4988235", "CT")
    lactose_tolerant = lct != "CC"
    
    # Celiac risk
    dq2 = genotypes.get("rs2187668", "CC")
    dq8 = genotypes.get("rs7775228", "TT")
    celiac_risk = "HIGH" if dq2 == "TT" else "MODERATE" if dq2 == "CT" or dq8 == "CC" else "LOW"
    
    # Histamine
    dao = genotypes.get("rs1049793", "CC")
    histamine_intolerant = dao == "GG"
    
    # Alcohol
    aldh2 = genotypes.get("rs671", "GG")
    adh1b = genotypes.get("rs1229984", "CC")
    alcohol_flush = aldh2 in ["AA", "GA"] or adh1b in ["TT", "CT"]
    
    # Bitter taste
    tas2r38 = genotypes.get("rs713598", "CG")
    bitter_taster = "super" if tas2r38 == "CC" else "non" if tas2r38 == "GG" else "medium"
    
    return {
        "apoe_genotype": apoe,
        "diet_recommendations": {
            "type": apoe_recs["diet_type"].value,
            "saturated_fat_tolerance": "LOW - strict limit" if sat_fat_sensitive else "Moderate",
            "carb_tolerance": carb_tolerance.value,
            "sodium_sensitivity": "HIGH - restrict sodium" if salt_sensitive else "Normal"
        },
        "food_intolerances": {
            "lactose": "Intolerant" if not lactose_tolerant else "Tolerant",
            "celiac_risk": celiac_risk,
            "histamine": "Likely intolerant" if histamine_intolerant else "Normal",
            "alcohol_flush": "Yes - flushing expected" if alcohol_flush else "Normal"
        },
        "taste_profile": {
            "bitter": bitter_taster,
            "cilantro_soapy": genotypes.get("rs72921001", "GG") == "AA"
        },
        "key_recommendations": apoe_recs["recommendations"],
        "markers_analyzed": sum(1 for rs in NUTRITION_COMPREHENSIVE_MARKERS if rs in genotypes),
    }

# Export
__all__ = [
    'NUTRITION_COMPREHENSIVE_MARKERS',
    'FAT_METABOLISM_MARKERS',
    'APOE_MARKERS',
    'APOE_DIET_RECOMMENDATIONS',
    'CARB_MARKERS',
    'SODIUM_MARKERS',
    'HISTAMINE_MARKERS',
    'OXALATE_MARKERS',
    'CELIAC_MARKERS',
    'LACTOSE_MARKERS',
    'FRUCTOSE_MARKERS',
    'ALCOHOL_MARKERS',
    'TASTE_MARKERS',
    'ToleranceLevel',
    'DietaryRecommendation',
    'determine_apoe_genotype',
    'generate_nutrition_report',
]
