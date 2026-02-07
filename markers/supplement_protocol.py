"""
Complete Supplement Protocol Generator v5.0

Genetic markers affecting nutrient metabolism and supplementation needs:
- Full methylation cycle
- Vitamin D pathway
- Vitamin A/B12/E metabolism
- Omega-3 conversion
- Iron metabolism
- Magnesium/Zinc status
- Antioxidant genes
- CoQ10 synthesis
- Choline metabolism

Each marker includes evidence-based supplementation recommendations.
"""

from enum import Enum
from typing import Dict, List, Any, Optional

class SupplementPriority(Enum):
    ESSENTIAL = "essential"      # Strong genetic need
    RECOMMENDED = "recommended"  # Moderate genetic support
    OPTIONAL = "optional"        # May benefit
    CAUTION = "caution"          # May need less or avoid

# =============================================================================
# FULL METHYLATION CYCLE
# =============================================================================

METHYLATION_MARKERS = {
    # MTHFR - Central folate enzyme
    "rs1801133": {
        "gene": "MTHFR",
        "variant": "C677T",
        "function": "5,10-MTHF to 5-MTHF conversion",
        "risk_allele": "T",
        "frequency": {"EUR": 0.35, "EAS": 0.35, "AFR": 0.10, "AMR": 0.50},
        "effect": {
            "TT": "~70% reduced enzyme activity",
            "CT": "~35% reduced enzyme activity",
            "CC": "Normal activity"
        },
        "category": "methylation",
        "evidence": "definitive",
        "pmid": ["8630491", "26647857", "16825332"],
        "supplement_impact": {
            "TT": {
                "methylfolate": SupplementPriority.ESSENTIAL,
                "folic_acid": SupplementPriority.CAUTION,
                "b12_methylcobalamin": SupplementPriority.RECOMMENDED,
                "riboflavin_b2": SupplementPriority.RECOMMENDED
            },
            "CT": {
                "methylfolate": SupplementPriority.RECOMMENDED,
                "folic_acid": SupplementPriority.OPTIONAL
            },
            "CC": {
                "folate": SupplementPriority.OPTIONAL
            }
        },
        "actionable": {
            "TT": [
                "Use methylfolate (L-5-MTHF) instead of folic acid",
                "Typical dose: 400-800mcg methylfolate daily",
                "Consider methylcobalamin (B12) 1000mcg",
                "Riboflavin (B2) 25-50mg may improve enzyme stability",
                "Check homocysteine level - treat if elevated",
                "Avoid folic acid fortified foods if possible",
                "Important during pregnancy - use prenatal with methylfolate"
            ],
            "CT": [
                "Consider methylfolate over folic acid",
                "Moderate dose: 400mcg methylfolate",
                "May use regular folic acid if methylfolate unavailable"
            ],
            "CC": [
                "No specific methylation support needed",
                "Standard folate intake adequate"
            ]
        }
    },
    "rs1801131": {
        "gene": "MTHFR",
        "variant": "A1298C",
        "function": "BH4 (tetrahydrobiopterin) synthesis pathway",
        "risk_allele": "C",
        "frequency": {"EUR": 0.33, "EAS": 0.20, "AFR": 0.15},
        "effect": {
            "CC": "~40% reduced enzyme activity for this pathway",
            "AC": "~20% reduced activity"
        },
        "category": "methylation",
        "evidence": "moderate",
        "pmid": ["10444342", "11684860"],
        "supplement_impact": {
            "CC": {
                "methylfolate": SupplementPriority.RECOMMENDED,
                "bh4_precursors": SupplementPriority.OPTIONAL
            }
        },
        "actionable": {
            "CC": [
                "Consider methylfolate supplementation",
                "Compound heterozygous (677CT + 1298AC) has additive effect",
                "BH4 pathway affects neurotransmitter synthesis"
            ]
        },
        "notes": "Affects BH4 rather than direct homocysteine metabolism"
    },
    
    # MTR - Methionine Synthase
    "rs1805087": {
        "gene": "MTR",
        "variant": "A2756G",
        "function": "B12-dependent homocysteine remethylation",
        "risk_allele": "G",
        "frequency": {"EUR": 0.20, "EAS": 0.15, "AFR": 0.35},
        "effect": {
            "GG": "Altered enzyme function - increased B12 requirements",
            "AG": "Mild effect"
        },
        "category": "methylation",
        "evidence": "moderate",
        "pmid": ["11114891", "21930180"],
        "supplement_impact": {
            "GG": {
                "b12_methylcobalamin": SupplementPriority.ESSENTIAL,
                "methylfolate": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "GG": [
                "Higher B12 requirements",
                "Use methylcobalamin or hydroxocobalamin forms",
                "Typical dose: 1000-2000mcg B12 daily",
                "Check B12 and MMA (methylmalonic acid) levels"
            ]
        }
    },
    
    # MTRR - Methionine Synthase Reductase
    "rs1801394": {
        "gene": "MTRR",
        "variant": "A66G",
        "function": "Regenerates MTR enzyme (B12 recycling)",
        "risk_allele": "G",
        "frequency": {"EUR": 0.50, "EAS": 0.25, "AFR": 0.30},
        "effect": {
            "GG": "Reduced B12 regeneration efficiency",
            "AG": "Mild effect"
        },
        "category": "methylation",
        "evidence": "moderate",
        "pmid": ["12518998", "17914124"],
        "supplement_impact": {
            "GG": {
                "b12_methylcobalamin": SupplementPriority.RECOMMENDED,
                "riboflavin_b2": SupplementPriority.OPTIONAL
            }
        },
        "actionable": {
            "GG": [
                "May have higher B12 requirements",
                "Consider active B12 forms (methylcobalamin)",
                "Combines with MTR variants for B12 needs"
            ]
        }
    },
    
    # BHMT - Betaine-Homocysteine Methyltransferase
    "rs3733890": {
        "gene": "BHMT",
        "variant": "R239Q",
        "function": "Alternative homocysteine remethylation (uses betaine)",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30, "EAS": 0.20},
        "effect": {
            "AA": "May affect betaine-dependent methylation"
        },
        "category": "methylation",
        "evidence": "moderate",
        "pmid": ["15531384", "16684769"],
        "supplement_impact": {
            "AA": {
                "tmg_betaine": SupplementPriority.OPTIONAL,
                "choline": SupplementPriority.OPTIONAL
            }
        },
        "actionable": {
            "AA": [
                "May benefit from TMG (betaine) supplementation",
                "Choline supplementation may help (converts to betaine)",
                "Alternative methylation pathway support"
            ]
        }
    },
    
    # CBS - Cystathionine Beta-Synthase
    "rs234706": {
        "gene": "CBS",
        "variant": "C699T",
        "function": "Transsulfuration pathway (homocysteine to cysteine)",
        "risk_allele": "T",
        "frequency": {"EUR": 0.45, "EAS": 0.55},
        "effect": {
            "TT": "Faster CBS activity - may deplete B12/folate faster"
        },
        "category": "methylation",
        "evidence": "moderate",
        "pmid": ["12777168", "17389630"],
        "supplement_impact": {
            "TT": {
                "b6_p5p": SupplementPriority.RECOMMENDED,
                "methylfolate": SupplementPriority.RECOMMENDED,
                "sulfur_caution": SupplementPriority.CAUTION
            }
        },
        "actionable": {
            "TT": [
                "CBS upregulation may increase sulfur/ammonia load",
                "May need MORE methylation support (B12, folate)",
                "B6 (P5P form) cofactor for CBS",
                "Some practitioners suggest limiting sulfur foods - evidence limited"
            ]
        },
        "notes": "Controversial - some claims not well supported"
    },
    
    # COMT - Catechol-O-Methyltransferase
    "rs4680": {
        "gene": "COMT",
        "variant": "Val158Met",
        "function": "Dopamine/catecholamine degradation",
        "risk_allele": "A",  # Met = slow
        "frequency": {"EUR": 0.48, "EAS": 0.30, "AFR": 0.40},
        "effect": {
            "AA": "Met/Met - Slow COMT (3-4x slower catecholamine breakdown)",
            "AG": "Val/Met - Intermediate",
            "GG": "Val/Val - Fast COMT"
        },
        "category": "methylation",
        "evidence": "definitive",
        "pmid": ["12716966", "17008817", "14517761"],
        "supplement_impact": {
            "AA": {
                "methylfolate": SupplementPriority.CAUTION,
                "sam_e": SupplementPriority.CAUTION,
                "magnesium": SupplementPriority.RECOMMENDED,
                "theanine": SupplementPriority.OPTIONAL
            },
            "GG": {
                "methylfolate": SupplementPriority.OPTIONAL,
                "sam_e": SupplementPriority.OPTIONAL,
                "tyrosine": SupplementPriority.OPTIONAL
            }
        },
        "actionable": {
            "AA": [
                "Slow COMT - higher dopamine/norepinephrine baseline",
                "May be MORE sensitive to methylation supplements",
                "Start methyl donors at LOW doses - can cause anxiety",
                "May benefit from magnesium for stress",
                "L-theanine may help balance catecholamines",
                "'Worrier' phenotype - detail-oriented but anxiety-prone"
            ],
            "GG": [
                "Fast COMT - lower baseline dopamine",
                "May tolerate higher doses of methyl donors",
                "'Warrior' phenotype - stress resilient",
                "May benefit from tyrosine for dopamine support"
            ],
            "AG": [
                "Intermediate COMT activity",
                "Balanced approach to methylation support"
            ]
        }
    },
}

# =============================================================================
# VITAMIN D PATHWAY
# =============================================================================

VITAMIN_D_MARKERS = {
    "rs2282679": {
        "gene": "GC",
        "variant": "Vitamin D binding protein",
        "function": "Transports vitamin D in blood",
        "risk_allele": "G",
        "frequency": {"EUR": 0.28, "AFR": 0.90, "EAS": 0.10},
        "effect": {
            "GG": "Lower vitamin D binding protein - lower 25(OH)D levels",
            "AG": "Intermediate effect"
        },
        "category": "vitamin_d",
        "evidence": "definitive",
        "pmid": ["20541252", "20418485", "23001564"],
        "supplement_impact": {
            "GG": {
                "vitamin_d3": SupplementPriority.ESSENTIAL,
                "vitamin_k2": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "GG": [
                "Higher vitamin D supplementation likely needed",
                "Target higher end of normal range (40-60 ng/mL)",
                "Typical dose: 2000-5000 IU D3 daily (test and adjust)",
                "Take vitamin K2 (MK-7) with D3 for bone/cardiovascular health",
                "Take with fat-containing meal for absorption"
            ]
        }
    },
    "rs12785878": {
        "gene": "DHCR7/NADSYN1",
        "variant": "Vitamin D synthesis",
        "function": "7-dehydrocholesterol reductase - D synthesis in skin",
        "risk_allele": "G",
        "frequency": {"EUR": 0.25, "AFR": 0.05},
        "effect": {
            "GG": "Reduced vitamin D synthesis from sunlight",
            "TG": "Intermediate"
        },
        "category": "vitamin_d",
        "evidence": "definitive",
        "pmid": ["20541252", "21129711"],
        "supplement_impact": {
            "GG": {
                "vitamin_d3": SupplementPriority.ESSENTIAL
            }
        },
        "actionable": {
            "GG": [
                "Reduced ability to synthesize vitamin D from sun",
                "More dependent on dietary/supplemental D",
                "Higher supplementation needs"
            ]
        }
    },
    "rs10741657": {
        "gene": "CYP2R1",
        "variant": "25-hydroxylase",
        "function": "Liver conversion of D to 25(OH)D",
        "risk_allele": "A",
        "frequency": {"EUR": 0.40, "AFR": 0.20},
        "effect": {
            "AA": "Reduced 25-hydroxylation of vitamin D"
        },
        "category": "vitamin_d",
        "evidence": "definitive",
        "pmid": ["20541252", "23001564"],
        "supplement_impact": {
            "AA": {
                "vitamin_d3": SupplementPriority.ESSENTIAL,
                "calcifediol": SupplementPriority.OPTIONAL
            }
        },
        "actionable": {
            "AA": [
                "May need higher D3 doses for same serum level",
                "Consider calcifediol (25-OH-D3) if not responding to D3",
                "Monitor 25(OH)D levels closely"
            ]
        }
    },
    "rs703842": {
        "gene": "CYP27B1",
        "variant": "1α-hydroxylase",
        "function": "Kidney activation to 1,25(OH)2D",
        "risk_allele": "T",
        "frequency": {"EUR": 0.35},
        "effect": "Affects active vitamin D production",
        "category": "vitamin_d",
        "evidence": "moderate",
        "pmid": ["20541252"]
    },
    # VDR variants
    "rs1544410": {
        "gene": "VDR",
        "variant": "BsmI (A>G)",
        "function": "Vitamin D receptor",
        "risk_allele": "G",
        "frequency": {"EUR": 0.40, "AFR": 0.30, "EAS": 0.10},
        "effect": "Affects VDR expression/function",
        "category": "vitamin_d",
        "evidence": "moderate",
        "pmid": ["12777168", "19343046"]
    },
    "rs731236": {
        "gene": "VDR",
        "variant": "TaqI (T>C)",
        "function": "Vitamin D receptor",
        "risk_allele": "C",
        "frequency": {"EUR": 0.40, "AFR": 0.30},
        "effect": "Affects VDR expression",
        "category": "vitamin_d",
        "evidence": "moderate",
        "pmid": ["12777168"]
    },
    "rs10735810": {
        "gene": "VDR",
        "variant": "FokI (C>T)",
        "function": "Vitamin D receptor start codon",
        "risk_allele": "T",
        "frequency": {"EUR": 0.35, "AFR": 0.25, "EAS": 0.50},
        "effect": {
            "TT": "Longer, less active VDR protein",
            "CC": "Shorter, more active VDR"
        },
        "category": "vitamin_d",
        "evidence": "strong",
        "pmid": ["19343046", "17968459"],
        "supplement_impact": {
            "TT": {
                "vitamin_d3": SupplementPriority.ESSENTIAL
            }
        },
        "actionable": {
            "TT": [
                "Less responsive to vitamin D",
                "May need higher levels for same effect",
                "Target higher end of optimal range"
            ]
        }
    },
}

# =============================================================================
# VITAMIN A (Beta-Carotene Conversion)
# =============================================================================

VITAMIN_A_MARKERS = {
    "rs7501331": {
        "gene": "BCMO1",
        "variant": "R267S",
        "function": "Beta-carotene to retinol conversion",
        "risk_allele": "T",
        "frequency": {"EUR": 0.24, "AFR": 0.35, "EAS": 0.20},
        "effect": {
            "TT": "32% reduced conversion efficiency",
            "CT": "~16% reduced"
        },
        "category": "vitamin_a",
        "evidence": "strong",
        "pmid": ["19103647", "22113863"],
        "supplement_impact": {
            "TT": {
                "retinol": SupplementPriority.ESSENTIAL,
                "beta_carotene": SupplementPriority.CAUTION
            }
        },
        "actionable": {
            "TT": [
                "Poor converter - beta-carotene is NOT adequate vitamin A source",
                "Need preformed retinol (liver, eggs, dairy, fish)",
                "If vegetarian: consider retinol supplement (retinyl palmitate)",
                "Carotenemia (orange skin) possible if high beta-carotene diet",
                "Critical for vegetarians/vegans with this genotype"
            ]
        }
    },
    "rs12934922": {
        "gene": "BCMO1",
        "variant": "A379V",
        "function": "Beta-carotene conversion",
        "risk_allele": "T",
        "frequency": {"EUR": 0.42, "AFR": 0.20},
        "effect": {
            "TT": "~48% reduced conversion when combined with R267S"
        },
        "category": "vitamin_a",
        "evidence": "strong",
        "pmid": ["19103647"],
        "supplement_impact": {
            "TT": {
                "retinol": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "TT": [
                "Additional reduction in beta-carotene conversion",
                "Combined with R267S can cause severe conversion deficit",
                "Preformed vitamin A important"
            ]
        }
    },
}

# =============================================================================
# VITAMIN B12
# =============================================================================

VITAMIN_B12_MARKERS = {
    "rs602662": {
        "gene": "FUT2",
        "variant": "Secretor status",
        "function": "B12 absorption in gut",
        "risk_allele": "A",
        "frequency": {"EUR": 0.45, "AFR": 0.25, "EAS": 0.45},
        "effect": {
            "AA": "Non-secretor - reduced B12 absorption, altered gut microbiome"
        },
        "category": "vitamin_b12",
        "evidence": "strong",
        "pmid": ["19197348", "21878437"],
        "supplement_impact": {
            "AA": {
                "b12_sublingual": SupplementPriority.ESSENTIAL,
                "probiotics": SupplementPriority.OPTIONAL
            }
        },
        "actionable": {
            "AA": [
                "Non-secretor status - lower B12 levels typical",
                "May benefit from sublingual B12 (bypasses absorption issues)",
                "Higher B12 doses may be needed",
                "Also affects gut microbiome and norovirus resistance"
            ]
        }
    },
    "rs1801198": {
        "gene": "TCN2",
        "variant": "P259R",
        "function": "Transcobalamin II - B12 cellular delivery",
        "risk_allele": "G",
        "frequency": {"EUR": 0.45, "AFR": 0.30, "EAS": 0.60},
        "effect": {
            "GG": "Reduced B12 delivery to cells despite normal serum levels"
        },
        "category": "vitamin_b12",
        "evidence": "strong",
        "pmid": ["11160079", "19578716"],
        "supplement_impact": {
            "GG": {
                "b12_methylcobalamin": SupplementPriority.ESSENTIAL
            }
        },
        "actionable": {
            "GG": [
                "B12 cellular delivery impaired",
                "May have normal serum B12 but functional deficiency",
                "Check MMA and homocysteine (better markers)",
                "May need higher B12 supplementation",
                "Sublingual/injectable may help"
            ]
        }
    },
}

# =============================================================================
# VITAMIN E
# =============================================================================

VITAMIN_E_MARKERS = {
    "rs6994076": {
        "gene": "TTPA",
        "variant": "Alpha-tocopherol transfer protein",
        "function": "Vitamin E transport and retention",
        "risk_allele": "T",
        "frequency": {"EUR": 0.30, "EAS": 0.20},
        "effect": {
            "TT": "Lower plasma vitamin E levels"
        },
        "category": "vitamin_e",
        "evidence": "moderate",
        "pmid": ["21878437", "23063622"],
        "supplement_impact": {
            "TT": {
                "vitamin_e": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "TT": [
                "May have lower vitamin E levels",
                "Consider mixed tocopherols supplement",
                "Vitamin E important for antioxidant defense"
            ]
        }
    },
}

# =============================================================================
# OMEGA-3 FATTY ACID METABOLISM
# =============================================================================

OMEGA3_MARKERS = {
    "rs174546": {
        "gene": "FADS1",
        "variant": "Fatty acid desaturase 1",
        "function": "ALA to EPA/DHA conversion",
        "risk_allele": "C",
        "frequency": {"EUR": 0.35, "AFR": 0.95, "EAS": 0.55},
        "effect": {
            "CC": "Efficient converter - high omega-6/low omega-3 levels",
            "TT": "Poor converter - lower omega-6 but needs preformed EPA/DHA"
        },
        "category": "omega3",
        "evidence": "definitive",
        "pmid": ["21829377", "26025071", "17056409"],
        "supplement_impact": {
            "TT": {
                "fish_oil": SupplementPriority.ESSENTIAL,
                "algal_oil": SupplementPriority.ESSENTIAL
            },
            "CC": {
                "fish_oil": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "TT": [
                "Poor ALA to EPA/DHA conversion (<5%)",
                "MUST get preformed EPA/DHA from fish or supplements",
                "Plant omega-3s (flax, chia) are insufficient",
                "Fish oil: 1000-2000mg EPA+DHA daily",
                "Vegetarians: algal oil for DHA"
            ],
            "CC": [
                "Efficient converter but produces more omega-6",
                "Still benefit from preformed EPA/DHA",
                "May want to limit omega-6 intake"
            ]
        }
    },
    "rs174547": {
        "gene": "FADS1",
        "variant": "FADS1 regulatory",
        "function": "Fatty acid desaturase activity",
        "risk_allele": "C",
        "frequency": {"EUR": 0.35, "AFR": 0.95},
        "effect": "Similar to rs174546",
        "category": "omega3",
        "evidence": "definitive",
        "pmid": ["21829377"]
    },
    "rs174537": {
        "gene": "FADS1",
        "variant": "FADS1/2 cluster",
        "function": "Affects both FADS1 and FADS2",
        "risk_allele": "G",
        "frequency": {"EUR": 0.35},
        "effect": "Part of conversion haplotype",
        "category": "omega3",
        "evidence": "definitive",
        "pmid": ["17056409"]
    },
    "rs1535": {
        "gene": "FADS2",
        "variant": "Delta-6 desaturase",
        "function": "First step in ALA conversion",
        "risk_allele": "G",
        "frequency": {"EUR": 0.35, "AFR": 0.95},
        "effect": "Affects omega-3 conversion efficiency",
        "category": "omega3",
        "evidence": "strong",
        "pmid": ["21829377"]
    },
    "rs3756963": {
        "gene": "ELOVL2",
        "variant": "Fatty acid elongase",
        "function": "EPA to DHA elongation",
        "risk_allele": "A",
        "frequency": {"EUR": 0.20},
        "effect": {
            "AA": "May have reduced EPA to DHA conversion"
        },
        "category": "omega3",
        "evidence": "moderate",
        "pmid": ["23362303"],
        "supplement_impact": {
            "AA": {
                "dha_specific": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "AA": [
                "May have difficulty converting EPA to DHA",
                "Ensure DHA specifically in supplement (not just EPA)",
                "DHA important for brain and eye health"
            ]
        }
    },
}

# =============================================================================
# IRON METABOLISM
# =============================================================================

IRON_MARKERS = {
    "rs1800562": {
        "gene": "HFE",
        "variant": "C282Y",
        "function": "Iron absorption regulation",
        "risk_allele": "A",
        "frequency": {"EUR": 0.06, "AFR": 0.001, "EAS": 0.001},
        "effect": {
            "AA": "Hereditary hemochromatosis - high iron overload risk",
            "AG": "Carrier - mild increased absorption"
        },
        "category": "iron",
        "evidence": "definitive",
        "pmid": ["8696333", "17100799"],
        "supplement_impact": {
            "AA": {
                "iron": SupplementPriority.CAUTION,
                "vitamin_c_with_meals": SupplementPriority.CAUTION
            },
            "AG": {
                "iron": SupplementPriority.CAUTION
            }
        },
        "actionable": {
            "AA": [
                "HIGH RISK hereditary hemochromatosis",
                "AVOID iron supplements",
                "AVOID vitamin C supplements with meals (increases iron absorption)",
                "Monitor ferritin and transferrin saturation regularly",
                "May need therapeutic phlebotomy",
                "Avoid cooking in cast iron",
                "Limit red meat intake"
            ],
            "AG": [
                "Carrier - monitor iron levels periodically",
                "Generally avoid iron supplements unless deficient"
            ]
        }
    },
    "rs1799945": {
        "gene": "HFE",
        "variant": "H63D",
        "function": "Iron absorption (milder)",
        "risk_allele": "G",
        "frequency": {"EUR": 0.14, "SAS": 0.10},
        "effect": {
            "GG": "Mild hemochromatosis risk",
            "AG": "Usually benign alone, risk if compound het with C282Y"
        },
        "category": "iron",
        "evidence": "definitive",
        "pmid": ["8696333", "17100799"],
        "supplement_impact": {
            "GG": {
                "iron": SupplementPriority.CAUTION
            }
        },
        "actionable": {
            "GG": [
                "Mild iron overload risk",
                "Monitor ferritin periodically"
            ],
            "compound_het": [
                "C282Y/H63D compound heterozygote - moderate risk",
                "Monitor iron studies annually"
            ]
        }
    },
    "rs855791": {
        "gene": "TMPRSS6",
        "variant": "A736V (matriptase-2)",
        "function": "Hepcidin regulation - iron absorption control",
        "risk_allele": "G",
        "frequency": {"EUR": 0.45, "EAS": 0.55, "AFR": 0.75},
        "effect": {
            "GG": "Lower iron levels, higher iron requirements"
        },
        "category": "iron",
        "evidence": "definitive",
        "pmid": ["19198612", "24019488"],
        "supplement_impact": {
            "GG": {
                "iron": SupplementPriority.OPTIONAL
            }
        },
        "actionable": {
            "GG": [
                "Genetically lower iron levels",
                "May need iron supplementation more often",
                "Monitor ferritin",
                "Important for menstruating women"
            ]
        }
    },
    "rs3811647": {
        "gene": "TF",
        "variant": "Transferrin variant",
        "function": "Iron transport protein",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30},
        "effect": "Affects iron transport capacity",
        "category": "iron",
        "evidence": "moderate",
        "pmid": ["19198612"]
    },
}

# =============================================================================
# MAGNESIUM
# =============================================================================

MAGNESIUM_MARKERS = {
    "rs11144134": {
        "gene": "TRPM6",
        "variant": "Magnesium channel",
        "function": "Intestinal and renal magnesium absorption",
        "risk_allele": "T",
        "frequency": {"EUR": 0.15},
        "effect": {
            "TT": "Reduced magnesium absorption"
        },
        "category": "magnesium",
        "evidence": "moderate",
        "pmid": ["20581029", "18628520"],
        "supplement_impact": {
            "TT": {
                "magnesium": SupplementPriority.ESSENTIAL
            }
        },
        "actionable": {
            "TT": [
                "Higher magnesium requirements",
                "Consider magnesium supplementation (glycinate, citrate)",
                "Typical dose: 200-400mg elemental magnesium",
                "Avoid magnesium oxide (poor absorption)",
                "Important for hundreds of enzymatic reactions"
            ]
        }
    },
    "rs7925578": {
        "gene": "CNNM2",
        "variant": "Magnesium transporter",
        "function": "Renal magnesium handling",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35},
        "effect": "Affects serum magnesium levels",
        "category": "magnesium",
        "evidence": "strong",
        "pmid": ["20581029"],
        "supplement_impact": {
            "AA": {
                "magnesium": SupplementPriority.RECOMMENDED
            }
        }
    },
}

# =============================================================================
# ZINC
# =============================================================================

ZINC_MARKERS = {
    "rs13266634": {
        "gene": "SLC30A8",
        "variant": "R325W",
        "function": "Zinc transporter in pancreatic beta cells",
        "risk_allele": "T",
        "frequency": {"EUR": 0.70, "EAS": 0.60, "AFR": 0.95},
        "effect": {
            "CC": "Higher T2D risk, altered zinc/insulin handling"
        },
        "category": "zinc",
        "evidence": "definitive",
        "pmid": ["17463249", "24136000"],
        "supplement_impact": {
            "CC": {
                "zinc": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "CC": [
                "Associated with T2D risk and zinc status",
                "May benefit from zinc supplementation",
                "Zinc important for insulin storage/release",
                "Typical dose: 15-30mg zinc daily"
            ]
        },
        "notes": "Paradoxically, W allele (T) is protective against T2D"
    },
}

# =============================================================================
# ANTIOXIDANT GENES
# =============================================================================

ANTIOXIDANT_MARKERS = {
    "rs4880": {
        "gene": "SOD2",
        "variant": "Ala16Val (MnSOD)",
        "function": "Mitochondrial superoxide dismutase",
        "risk_allele": "T",  # Val allele
        "frequency": {"EUR": 0.47, "EAS": 0.10, "AFR": 0.45},
        "effect": {
            "TT": "Val/Val - enzyme targeted to mitochondria inefficiently",
            "CC": "Ala/Ala - efficient targeting but may produce excess H2O2"
        },
        "category": "antioxidant",
        "evidence": "strong",
        "pmid": ["12618587", "15534163"],
        "supplement_impact": {
            "TT": {
                "coq10": SupplementPriority.RECOMMENDED,
                "mito_support": SupplementPriority.RECOMMENDED
            },
            "CC": {
                "catalase_support": SupplementPriority.OPTIONAL
            }
        },
        "actionable": {
            "TT": [
                "Reduced mitochondrial antioxidant defense",
                "CoQ10 supplementation may help (100-200mg)",
                "Avoid excessive iron (oxidative stress)",
                "Mitochondrial support important"
            ],
            "CC": [
                "Efficient SOD2 but may produce excess H2O2",
                "Ensure adequate catalase/GPX cofactors (selenium)",
                "Balance antioxidant intake"
            ]
        }
    },
    "rs1001179": {
        "gene": "CAT",
        "variant": "Catalase promoter",
        "function": "Hydrogen peroxide breakdown",
        "risk_allele": "T",
        "frequency": {"EUR": 0.20, "AFR": 0.10},
        "effect": {
            "TT": "Lower catalase expression"
        },
        "category": "antioxidant",
        "evidence": "moderate",
        "pmid": ["10749980"],
        "supplement_impact": {
            "TT": {
                "glutathione_support": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "TT": [
                "Lower catalase activity",
                "Support glutathione peroxidase pathway as backup",
                "Ensure adequate selenium intake",
                "NAC may help support antioxidant systems"
            ]
        }
    },
    "rs1050450": {
        "gene": "GPX1",
        "variant": "Pro198Leu",
        "function": "Glutathione peroxidase - selenium enzyme",
        "risk_allele": "T",
        "frequency": {"EUR": 0.28, "EAS": 0.20, "AFR": 0.40},
        "effect": {
            "TT": "Reduced GPX1 activity, lower selenium utilization"
        },
        "category": "antioxidant",
        "evidence": "moderate",
        "pmid": ["12072403", "18174241"],
        "supplement_impact": {
            "TT": {
                "selenium": SupplementPriority.ESSENTIAL,
                "nac": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "TT": [
                "Reduced glutathione peroxidase activity",
                "Selenium supplementation important (200mcg)",
                "NAC (600-1200mg) supports glutathione",
                "Critical antioxidant enzyme"
            ]
        }
    },
    "rs1800566": {
        "gene": "NQO1",
        "variant": "Pro187Ser",
        "function": "NAD(P)H quinone oxidoreductase - detox and CoQ10 reduction",
        "risk_allele": "T",
        "frequency": {"EUR": 0.20, "EAS": 0.45, "AFR": 0.04},
        "effect": {
            "TT": "Complete loss of NQO1 activity",
            "CT": "~50% reduced activity"
        },
        "category": "antioxidant",
        "evidence": "strong",
        "pmid": ["12692552", "16249184"],
        "supplement_impact": {
            "TT": {
                "ubiquinol": SupplementPriority.ESSENTIAL,
                "nac": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "TT": [
                "Cannot reduce CoQ10 to ubiquinol efficiently",
                "MUST use ubiquinol form (not ubiquinone) for CoQ10",
                "NAC may help compensate",
                "Higher cancer risk - avoid benzene exposure"
            ],
            "CT": [
                "Prefer ubiquinol over ubiquinone for CoQ10"
            ]
        }
    },
}

# =============================================================================
# CoQ10
# =============================================================================

COQ10_MARKERS = {
    "rs1800566": {
        "gene": "NQO1",
        "variant": "Pro187Ser",
        "reference": "See ANTIOXIDANT_MARKERS"
    },
    "rs9380585": {
        "gene": "PDSS2",
        "variant": "CoQ10 synthesis",
        "function": "Prenyl diphosphate synthase - CoQ10 biosynthesis",
        "risk_allele": "A",
        "frequency": {"EUR": 0.25},
        "effect": {
            "AA": "May have reduced CoQ10 synthesis"
        },
        "category": "coq10",
        "evidence": "moderate",
        "pmid": ["21878437"],
        "supplement_impact": {
            "AA": {
                "coq10": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "AA": [
                "May have lower endogenous CoQ10 levels",
                "Consider CoQ10 supplementation (100-200mg)",
                "Ubiquinol form preferred if NQO1 TT"
            ]
        }
    },
}

# =============================================================================
# CHOLINE
# =============================================================================

CHOLINE_MARKERS = {
    "rs12325817": {
        "gene": "PEMT",
        "variant": "G5465A",
        "function": "Phosphatidylethanolamine N-methyltransferase - endogenous choline synthesis",
        "risk_allele": "C",
        "frequency": {"EUR": 0.75, "EAS": 0.80, "AFR": 0.65},
        "effect": {
            "CC": "Reduced ability to synthesize choline from SAMe/folate",
            "particularly_important": "In women (estrogen normally upregulates PEMT)"
        },
        "category": "choline",
        "evidence": "strong",
        "pmid": ["16702265", "22045291"],
        "supplement_impact": {
            "CC": {
                "choline": SupplementPriority.ESSENTIAL,
                "eggs": SupplementPriority.RECOMMENDED
            }
        },
        "actionable": {
            "CC": [
                "Higher dietary choline requirements",
                "Women especially affected (estrogen normally compensates)",
                "Eggs are best dietary source (~150mg/egg)",
                "Consider choline supplement (500mg)",
                "Important during pregnancy for fetal brain development",
                "May develop fatty liver on low-choline diet"
            ]
        }
    },
}

# =============================================================================
# COMBINE ALL SUPPLEMENT-RELATED MARKERS
# =============================================================================

SUPPLEMENT_MARKERS = {
    **METHYLATION_MARKERS,
    **VITAMIN_D_MARKERS,
    **VITAMIN_A_MARKERS,
    **VITAMIN_B12_MARKERS,
    **VITAMIN_E_MARKERS,
    **OMEGA3_MARKERS,
    **IRON_MARKERS,
    **MAGNESIUM_MARKERS,
    **ZINC_MARKERS,
    **ANTIOXIDANT_MARKERS,
    **COQ10_MARKERS,
    **CHOLINE_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def generate_supplement_protocol(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate personalized supplement protocol based on genetics."""
    
    recommendations = {
        "essential": [],
        "recommended": [],
        "optional": [],
        "caution": [],
        "summary": {}
    }
    
    # Methylation Analysis
    mthfr_677 = genotypes.get("rs1801133", "CC")
    mthfr_1298 = genotypes.get("rs1801131", "AA")
    comt = genotypes.get("rs4680", "AG")
    
    if mthfr_677 == "TT":
        recommendations["essential"].append({
            "supplement": "Methylfolate (L-5-MTHF)",
            "dose": "400-800mcg daily",
            "reason": "MTHFR C677T TT - significantly reduced enzyme activity",
            "pmid": "8630491"
        })
        recommendations["recommended"].append({
            "supplement": "Methylcobalamin (B12)",
            "dose": "1000mcg daily",
            "reason": "Supports methylation cycle with MTHFR variant"
        })
        recommendations["caution"].append({
            "supplement": "Folic acid",
            "reason": "May not convert well; use methylfolate instead"
        })
    elif mthfr_677 == "CT":
        recommendations["recommended"].append({
            "supplement": "Methylfolate",
            "dose": "400mcg daily",
            "reason": "MTHFR C677T CT - moderately reduced activity"
        })
    
    # COMT considerations
    if comt == "AA":  # Met/Met
        recommendations["caution"].append({
            "supplement": "High-dose methylation supplements",
            "reason": "Slow COMT - may increase anxiety. Start low, go slow."
        })
        recommendations["optional"].append({
            "supplement": "Magnesium (glycinate)",
            "dose": "200-400mg",
            "reason": "May help with stress/anxiety in slow COMT"
        })
    
    # Vitamin D
    gc = genotypes.get("rs2282679", "AA")
    vdr_fok = genotypes.get("rs10735810", "CC")
    
    if gc == "GG" or vdr_fok == "TT":
        recommendations["essential"].append({
            "supplement": "Vitamin D3",
            "dose": "2000-5000 IU daily (test and adjust)",
            "reason": f"{'GC GG: lower DBP, ' if gc == 'GG' else ''}{'VDR FokI TT: less responsive VDR' if vdr_fok == 'TT' else ''}",
            "pmid": "20541252"
        })
        recommendations["recommended"].append({
            "supplement": "Vitamin K2 (MK-7)",
            "dose": "100-200mcg",
            "reason": "Synergistic with D3 for bone and cardiovascular health"
        })
    
    # Beta-carotene conversion
    bcmo1 = genotypes.get("rs7501331", "CC")
    if bcmo1 == "TT":
        recommendations["essential"].append({
            "supplement": "Preformed Vitamin A (retinol)",
            "dose": "2500-5000 IU",
            "reason": "BCMO1 TT - poor beta-carotene converter",
            "pmid": "19103647"
        })
        recommendations["caution"].append({
            "supplement": "Beta-carotene as sole vitamin A source",
            "reason": "Inefficient conversion - need preformed retinol"
        })
    
    # B12
    fut2 = genotypes.get("rs602662", "GG")
    tcn2 = genotypes.get("rs1801198", "CC")
    
    if fut2 == "AA" or tcn2 == "GG":
        recommendations["essential"].append({
            "supplement": "Methylcobalamin (B12) sublingual",
            "dose": "1000-2000mcg",
            "reason": f"{'FUT2 non-secretor: absorption issues, ' if fut2 == 'AA' else ''}{'TCN2 GG: cellular delivery impaired' if tcn2 == 'GG' else ''}"
        })
    
    # Omega-3
    fads1 = genotypes.get("rs174546", "CT")
    if fads1 == "TT":
        recommendations["essential"].append({
            "supplement": "Fish oil or Algal oil",
            "dose": "1000-2000mg EPA+DHA",
            "reason": "FADS1 TT - poor ALA to EPA/DHA converter",
            "pmid": "21829377"
        })
    
    # Iron - Hemochromatosis
    hfe_c282y = genotypes.get("rs1800562", "GG")
    hfe_h63d = genotypes.get("rs1799945", "CC")
    
    if hfe_c282y == "AA":
        recommendations["caution"].append({
            "supplement": "Iron supplements",
            "reason": "HFE C282Y homozygous - hemochromatosis risk. AVOID iron."
        })
        recommendations["caution"].append({
            "supplement": "Vitamin C with meals",
            "reason": "Increases iron absorption - avoid with hemochromatosis risk"
        })
    elif hfe_c282y == "AG" or hfe_h63d == "GG":
        recommendations["caution"].append({
            "supplement": "Iron supplements",
            "reason": "HFE variant carrier - monitor iron levels, avoid supplements unless deficient"
        })
    
    # Antioxidants
    nqo1 = genotypes.get("rs1800566", "CC")
    sod2 = genotypes.get("rs4880", "CC")
    gpx1 = genotypes.get("rs1050450", "CC")
    
    if nqo1 == "TT":
        recommendations["essential"].append({
            "supplement": "Ubiquinol (reduced CoQ10)",
            "dose": "100-200mg",
            "reason": "NQO1 TT - cannot reduce CoQ10. Must use ubiquinol form.",
            "pmid": "12692552"
        })
    elif nqo1 == "CT":
        recommendations["recommended"].append({
            "supplement": "Ubiquinol (preferred over ubiquinone)",
            "dose": "100mg",
            "reason": "NQO1 CT - reduced CoQ10 conversion ability"
        })
    
    if gpx1 == "TT":
        recommendations["essential"].append({
            "supplement": "Selenium",
            "dose": "200mcg",
            "reason": "GPX1 TT - reduced glutathione peroxidase activity",
            "pmid": "12072403"
        })
        recommendations["recommended"].append({
            "supplement": "NAC (N-acetyl cysteine)",
            "dose": "600-1200mg",
            "reason": "Supports glutathione with GPX1 variant"
        })
    
    # Choline
    pemt = genotypes.get("rs12325817", "GG")
    if pemt == "CC":
        recommendations["essential"].append({
            "supplement": "Choline (or eggs)",
            "dose": "500mg or 3+ eggs daily",
            "reason": "PEMT CC - reduced endogenous choline synthesis",
            "pmid": "16702265"
        })
    
    # Magnesium
    trpm6 = genotypes.get("rs11144134", "CC")
    if trpm6 == "TT":
        recommendations["essential"].append({
            "supplement": "Magnesium (glycinate/citrate)",
            "dose": "200-400mg",
            "reason": "TRPM6 TT - reduced magnesium absorption"
        })
    
    # Generate summary
    recommendations["summary"] = {
        "methylation_status": "Impaired" if mthfr_677 == "TT" else "Moderate" if mthfr_677 == "CT" else "Normal",
        "comt_type": "Slow (worrier)" if comt == "AA" else "Fast (warrior)" if comt == "GG" else "Intermediate",
        "vitamin_d_needs": "High" if gc == "GG" or vdr_fok == "TT" else "Normal",
        "beta_carotene_converter": "Poor" if bcmo1 == "TT" else "Normal",
        "omega3_converter": "Poor" if fads1 == "TT" else "Normal",
        "iron_status": "Hemochromatosis risk" if hfe_c282y in ["AA", "AG"] else "Normal",
        "coq10_form_needed": "Ubiquinol required" if nqo1 == "TT" else "Either form OK",
        "total_essential_supplements": len(recommendations["essential"]),
        "total_recommended_supplements": len(recommendations["recommended"]),
        "warnings": len(recommendations["caution"])
    }
    
    return recommendations

# Export
__all__ = [
    'SUPPLEMENT_MARKERS',
    'METHYLATION_MARKERS',
    'VITAMIN_D_MARKERS',
    'VITAMIN_A_MARKERS',
    'VITAMIN_B12_MARKERS',
    'VITAMIN_E_MARKERS',
    'OMEGA3_MARKERS',
    'IRON_MARKERS',
    'MAGNESIUM_MARKERS',
    'ZINC_MARKERS',
    'ANTIOXIDANT_MARKERS',
    'COQ10_MARKERS',
    'CHOLINE_MARKERS',
    'SupplementPriority',
    'generate_supplement_protocol',
]
