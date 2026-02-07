"""
Comprehensive Athletic Genetics Panel v5.0

Complete coverage of:
- Power vs Endurance genetics
- VO2max potential
- Muscle fiber type
- Lactate metabolism
- Recovery genetics
- Injury risk (tendon, ligament, bone)
- ACL and Achilles risk
- Creatine response
- Caffeine ergogenic effect
- Heat tolerance
- Altitude adaptation
- Concussion recovery
- Hydration genetics

All markers with PMID references and sport-specific recommendations.
"""

from enum import Enum
from typing import Dict, List, Any, Optional

class AthleticProfile(Enum):
    POWER_DOMINANT = "power_dominant"        # Sprinter, weightlifter
    POWER_ENDURANCE = "power_endurance"      # Mixed profile
    BALANCED = "balanced"                    # Versatile
    ENDURANCE_POWER = "endurance_power"      # Mixed with endurance lean
    ENDURANCE_DOMINANT = "endurance_dominant" # Marathon, triathlon

class RecoverySpeed(Enum):
    FAST = "fast"           # <24h to baseline
    NORMAL = "normal"       # 24-48h
    SLOW = "slow"           # 48-72h+

class InjuryRiskLevel(Enum):
    LOW = "low"
    MODERATE = "moderate"
    ELEVATED = "elevated"
    HIGH = "high"

# =============================================================================
# POWER / MUSCLE FIBER TYPE
# =============================================================================

POWER_MARKERS = {
    "rs1815739": {
        "gene": "ACTN3",
        "variant": "R577X",
        "function": "Alpha-actinin-3 - fast-twitch muscle fiber protein",
        "risk_allele": "T",  # X allele = stop codon
        "frequency": {"EUR": 0.44, "EAS": 0.55, "AFR": 0.15},
        "effect": {
            "CC": "RR - Full ACTN3 expression. Power/sprint advantage. ~18% of population.",
            "CT": "RX - Intermediate. Good for power sports.",
            "TT": "XX - No ACTN3. Endurance advantage. Elite sprinters virtually never XX."
        },
        "category": "power",
        "evidence": "definitive",
        "pmid": ["12879365", "18043716", "20308985"],
        "sport_implications": {
            "CC": ["Sprint", "Powerlifting", "Shot put", "Football (lineman)", "Weightlifting"],
            "CT": ["Soccer", "Basketball", "Swimming", "Martial arts", "All-round sports"],
            "TT": ["Marathon", "Triathlon", "Cycling", "Cross-country", "Ultrarunning"]
        },
        "actionable": {
            "CC": [
                "Power/sprint genetic advantage",
                "Higher fast-twitch fiber proportion",
                "Responds well to strength/power training",
                "May not excel at ultra-endurance naturally",
                "Consider sports requiring explosive power"
            ],
            "TT": [
                "Endurance genetic profile",
                "Higher slow-twitch fiber proportion",
                "Natural endurance ability",
                "May need more power-focused training if doing sprint sports",
                "Elite sprinters are never XX genotype"
            ]
        }
    },
    "rs17602729": {
        "gene": "AMPD1",
        "variant": "Q12X (Gln12Ter)",
        "function": "AMP deaminase - energy metabolism in muscle",
        "risk_allele": "A",
        "frequency": {"EUR": 0.12, "AFR": 0.02},
        "effect": {
            "AA": "AMPD deficiency - exercise intolerance, cramps",
            "AG": "Carrier - may have mild exercise intolerance",
            "GG": "Normal AMPD function"
        },
        "category": "power",
        "evidence": "strong",
        "pmid": ["12403786", "18997033"],
        "actionable": {
            "AA": [
                "May experience exercise-induced muscle cramps/pain",
                "Longer warm-up recommended",
                "May have reduced high-intensity capacity",
                "Not incompatible with exercise, just may need adaptation"
            ]
        }
    },
}

# =============================================================================
# ENDURANCE / VO2MAX
# =============================================================================

ENDURANCE_MARKERS = {
    "rs4340": {
        "gene": "ACE",
        "variant": "Insertion/Deletion (I/D)",
        "function": "Angiotensin-converting enzyme levels",
        "risk_allele": "I",  # I allele = endurance
        "frequency": {"EUR": 0.48, "EAS": 0.60, "AFR": 0.40},
        "effect": {
            "II": "Low ACE levels. Endurance advantage. Elite in mountaineering, distance running.",
            "ID": "Intermediate ACE. Balanced performance.",
            "DD": "High ACE levels. Power/strength advantage. Blood pressure sensitivity."
        },
        "category": "endurance",
        "evidence": "strong",
        "pmid": ["8675673", "10917533", "16565348"],
        "sport_implications": {
            "II": ["Marathon", "Mountaineering", "Rowing", "Cycling", "Ultra-endurance"],
            "ID": ["Swimming", "Soccer", "Tennis", "Triathlons"],
            "DD": ["Powerlifting", "Sprinting", "Short-distance events"]
        },
        "actionable": {
            "II": [
                "Endurance genetic advantage",
                "Efficient oxygen utilization",
                "Better altitude adaptation",
                "Natural long-distance ability",
                "May need dedicated strength work for power sports"
            ],
            "DD": [
                "Power/strength genetic advantage",
                "May excel at short, intense efforts",
                "Salt-sensitive blood pressure (monitor if hypertensive)",
                "May need dedicated endurance training for distance events"
            ]
        }
    },
    "rs8192678": {
        "gene": "PPARGC1A",
        "variant": "Gly482Ser (PGC-1α)",
        "function": "Master regulator of mitochondrial biogenesis",
        "risk_allele": "A",  # Ser = reduced function
        "frequency": {"EUR": 0.35, "EAS": 0.40, "AFR": 0.15},
        "effect": {
            "GG": "Gly/Gly - Optimal PGC-1α function. Higher VO2max potential.",
            "GA": "Intermediate mitochondrial response",
            "AA": "Ser/Ser - Reduced mitochondrial biogenesis response to training"
        },
        "category": "endurance",
        "evidence": "strong",
        "pmid": ["14586024", "15677319", "18500963"],
        "actionable": {
            "GG": [
                "High VO2max trainability",
                "Responds well to endurance training",
                "Better mitochondrial adaptation to exercise"
            ],
            "AA": [
                "May need more training stimulus for VO2max gains",
                "Higher intensity intervals may be more effective",
                "Consider longer training blocks"
            ]
        }
    },
    "rs2010963": {
        "gene": "VEGFA",
        "variant": "VEGF -634G>C",
        "function": "Vascular endothelial growth factor - angiogenesis",
        "risk_allele": "C",
        "frequency": {"EUR": 0.30, "EAS": 0.40},
        "effect": {
            "GG": "Higher VEGF levels - better capillarization",
            "CC": "Lower VEGF - reduced angiogenic response"
        },
        "category": "endurance",
        "evidence": "moderate",
        "pmid": ["16314865", "15772360"],
        "actionable": {
            "GG": [
                "Better capillary response to training",
                "Good endurance adaptation capacity"
            ]
        }
    },
    "rs11549465": {
        "gene": "HIF1A",
        "variant": "Pro582Ser",
        "function": "Hypoxia-inducible factor - altitude/hypoxia response",
        "risk_allele": "T",  # Ser allele
        "frequency": {"EUR": 0.05, "EAS": 0.02},
        "effect": {
            "CT/TT": "Enhanced hypoxia response - potential altitude advantage"
        },
        "category": "endurance",
        "evidence": "moderate",
        "pmid": ["14507923", "17000709"],
        "actionable": {
            "CT": [
                "May have enhanced altitude adaptation",
                "Could benefit from altitude training"
            ]
        }
    },
    "rs12594956": {
        "gene": "NRF1",
        "variant": "Nuclear respiratory factor 1",
        "function": "Mitochondrial gene expression",
        "risk_allele": "A",
        "frequency": {"EUR": 0.40},
        "effect": "Affects mitochondrial biogenesis",
        "category": "endurance",
        "evidence": "moderate",
        "pmid": ["18500963"]
    },
    "rs8192675": {
        "gene": "CKM",
        "variant": "Creatine kinase, muscle",
        "function": "Energy transfer in muscle",
        "risk_allele": "T",
        "frequency": {"EUR": 0.30},
        "effect": "Associated with VO2max trainability",
        "category": "endurance",
        "evidence": "moderate",
        "pmid": ["22971578"],
        "actionable": {
            "TT": [
                "Higher VO2max training response",
                "Good aerobic trainability"
            ]
        }
    },
}

# =============================================================================
# LACTATE METABOLISM
# =============================================================================

LACTATE_MARKERS = {
    "rs1049434": {
        "gene": "SLC16A1",
        "variant": "MCT1 (A1470T)",
        "function": "Monocarboxylate transporter 1 - lactate shuttling",
        "risk_allele": "T",
        "frequency": {"EUR": 0.40, "EAS": 0.30, "AFR": 0.20},
        "effect": {
            "AA": "Efficient lactate clearance - higher lactate threshold",
            "TT": "Reduced lactate transport - earlier fatigue at threshold",
            "AT": "Intermediate lactate handling"
        },
        "category": "lactate",
        "evidence": "strong",
        "pmid": ["16009822", "19574090"],
        "actionable": {
            "AA": [
                "Efficient lactate clearance",
                "Higher lactate threshold potential",
                "Can sustain higher intensities longer",
                "Responds well to threshold training"
            ],
            "TT": [
                "May fatigue earlier at high intensities",
                "Focus on building lactate tolerance",
                "Threshold training especially important"
            ]
        }
    },
    "rs3849364": {
        "gene": "LACTB",
        "variant": "Lactamase beta",
        "function": "Mitochondrial protein affecting lactate metabolism",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30},
        "effect": "Associated with endurance performance",
        "category": "lactate",
        "evidence": "moderate",
        "pmid": ["19997062"]
    },
}

# =============================================================================
# RECOVERY
# =============================================================================

RECOVERY_MARKERS = {
    "rs1800795": {
        "gene": "IL6",
        "variant": "IL-6 -174G>C",
        "function": "Interleukin-6 - inflammatory response",
        "risk_allele": "C",
        "frequency": {"EUR": 0.45, "EAS": 0.05, "AFR": 0.10},
        "effect": {
            "CC": "Higher IL-6 response - stronger but prolonged inflammation",
            "GG": "Lower IL-6 - faster recovery from exercise",
            "GC": "Intermediate inflammatory response"
        },
        "category": "recovery",
        "evidence": "strong",
        "pmid": ["15616363", "15781034", "23632419"],
        "recovery_effect": {
            "CC": RecoverySpeed.SLOW,
            "GC": RecoverySpeed.NORMAL,
            "GG": RecoverySpeed.FAST
        },
        "actionable": {
            "CC": [
                "Stronger inflammatory response to training",
                "May need longer recovery between hard sessions",
                "Anti-inflammatory nutrition important (omega-3, tart cherry)",
                "Sleep especially critical for recovery",
                "Consider 48-72h between intense sessions"
            ],
            "GG": [
                "Efficient recovery from training",
                "Can tolerate higher training frequency",
                "May need higher training stimulus for adaptation"
            ]
        }
    },
    "rs1800629": {
        "gene": "TNF",
        "variant": "TNF-α -308G>A",
        "function": "Tumor necrosis factor alpha - inflammation",
        "risk_allele": "A",
        "frequency": {"EUR": 0.15, "AFR": 0.10},
        "effect": {
            "AA": "Very high TNF-α - prolonged inflammation",
            "GA": "Elevated inflammatory response",
            "GG": "Normal inflammation"
        },
        "category": "recovery",
        "evidence": "strong",
        "pmid": ["15616363", "12904686"],
        "recovery_effect": {
            "AA": RecoverySpeed.SLOW,
            "GA": RecoverySpeed.SLOW,
            "GG": RecoverySpeed.NORMAL
        },
        "actionable": {
            "AA": [
                "High inflammatory phenotype",
                "Extended recovery periods essential",
                "Anti-inflammatory strategies critical"
            ]
        }
    },
    "rs1205": {
        "gene": "CRP",
        "variant": "C-reactive protein variant",
        "function": "Baseline inflammation marker",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35, "AFR": 0.40},
        "effect": {
            "AA": "Higher baseline CRP - more inflammation-prone"
        },
        "category": "recovery",
        "evidence": "strong",
        "pmid": ["17157862", "19270852"]
    },
}

# =============================================================================
# INJURY RISK - TENDONS/LIGAMENTS
# =============================================================================

TENDON_MARKERS = {
    "rs1800012": {
        "gene": "COL1A1",
        "variant": "Sp1 binding site (G>T)",
        "function": "Collagen type 1 alpha 1 - structural protein",
        "risk_allele": "T",
        "frequency": {"EUR": 0.20, "AFR": 0.10, "EAS": 0.05},
        "effect": {
            "TT": "Reduced collagen strength - HIGHER injury risk",
            "GT": "Intermediate - moderate injury risk",
            "GG": "Normal collagen - lower injury risk"
        },
        "category": "injury",
        "evidence": "strong",
        "pmid": ["16477526", "18442638", "20460586"],
        "injury_risk": {
            "TT": InjuryRiskLevel.HIGH,
            "GT": InjuryRiskLevel.ELEVATED,
            "GG": InjuryRiskLevel.LOW
        },
        "injuries_affected": ["ACL rupture", "Achilles tendinopathy", "Rotator cuff", "General tendon injuries"],
        "actionable": {
            "TT": [
                "HIGHER soft tissue injury risk",
                "Extended warm-up essential",
                "Progressive loading critical",
                "May benefit from collagen supplementation (10-15g + vitamin C)",
                "Consider injury prevention exercises",
                "Don't ignore minor pain signals"
            ],
            "GT": [
                "Moderate injury risk",
                "Good warm-up and progression important"
            ]
        }
    },
    "rs12722": {
        "gene": "COL5A1",
        "variant": "C/T (3' UTR)",
        "function": "Collagen type V alpha 1 - regulates fibril assembly",
        "risk_allele": "T",
        "frequency": {"EUR": 0.45, "EAS": 0.35, "AFR": 0.30},
        "effect": {
            "TT": "Associated with ACL rupture risk",
            "CC": "Protective for ACL, but may be less flexible"
        },
        "category": "injury",
        "evidence": "strong",
        "pmid": ["15689377", "20460586", "23054209"],
        "injury_risk": {
            "TT": InjuryRiskLevel.ELEVATED,
            "CT": InjuryRiskLevel.MODERATE,
            "CC": InjuryRiskLevel.LOW
        },
        "injuries_affected": ["ACL rupture", "Achilles tendinopathy", "Chronic Achilles problems"],
        "actionable": {
            "TT": [
                "Elevated ACL and Achilles injury risk",
                "Neuromuscular training programs recommended",
                "Landing mechanics training important",
                "Consider proprioception exercises"
            ]
        }
    },
    "rs240736": {
        "gene": "COL12A1",
        "variant": "Collagen XII variant",
        "function": "Collagen type XII - fibril organization",
        "risk_allele": "A",
        "frequency": {"EUR": 0.20},
        "effect": "Associated with ACL rupture risk",
        "category": "injury",
        "evidence": "moderate",
        "pmid": ["20460586", "22525324"],
        "injury_risk": {
            "AA": InjuryRiskLevel.ELEVATED
        }
    },
    "rs13321": {
        "gene": "TNC",
        "variant": "Tenascin-C GT repeat proxy",
        "function": "Tenascin-C - ECM protein in tendons",
        "risk_allele": "T",
        "frequency": {"EUR": 0.25},
        "effect": {
            "TT": "Associated with Achilles tendinopathy"
        },
        "category": "injury",
        "evidence": "strong",
        "pmid": ["15689377", "18927258"],
        "injuries_affected": ["Achilles tendinopathy", "Chronic Achilles problems"],
        "actionable": {
            "TT": [
                "Higher Achilles injury risk",
                "Careful progression with running volume",
                "Eccentric heel drops for prevention"
            ]
        }
    },
    "rs679620": {
        "gene": "MMP3",
        "variant": "5A/6A promoter",
        "function": "Matrix metalloproteinase 3 - tissue remodeling",
        "risk_allele": "T",  # 5A allele
        "frequency": {"EUR": 0.50, "EAS": 0.20},
        "effect": {
            "TT": "5A/5A - Higher MMP3, increased tendon breakdown",
            "CC": "6A/6A - Lower MMP3, may have stiffer tendons"
        },
        "category": "injury",
        "evidence": "moderate",
        "pmid": ["18927258", "19815014"],
        "injuries_affected": ["Achilles tendinopathy"]
    },
}

# =============================================================================
# BONE DENSITY
# =============================================================================

BONE_MARKERS = {
    "rs2234693": {
        "gene": "ESR1",
        "variant": "Estrogen receptor alpha PvuII",
        "function": "Estrogen receptor - bone density regulation",
        "risk_allele": "C",
        "frequency": {"EUR": 0.45, "EAS": 0.50, "AFR": 0.30},
        "effect": {
            "CC": "May have lower bone density response to exercise",
            "TT": "Better bone density adaptation"
        },
        "category": "bone",
        "evidence": "strong",
        "pmid": ["15014266", "15240626"],
        "actionable": {
            "CC": [
                "May need more impact exercise for bone density",
                "Ensure adequate calcium and vitamin D",
                "Weight-bearing exercise important"
            ]
        }
    },
    "rs1800012": {
        "gene": "COL1A1",
        "reference": "See TENDON_MARKERS",
        "effect": "Also affects bone density"
    },
    "rs1544410": {
        "gene": "VDR",
        "variant": "BsmI",
        "function": "Vitamin D receptor - calcium absorption, bone metabolism",
        "risk_allele": "G",
        "frequency": {"EUR": 0.40, "AFR": 0.30, "EAS": 0.10},
        "effect": "Affects bone mineral density",
        "category": "bone",
        "evidence": "moderate",
        "pmid": ["12777168", "14581593"],
        "actionable": {
            "GG": [
                "May have lower bone density",
                "Vitamin D optimization important",
                "Weight-bearing exercise recommended"
            ]
        }
    },
}

# =============================================================================
# CREATINE RESPONSE
# =============================================================================

CREATINE_MARKERS = {
    "rs1815739": {
        "gene": "ACTN3",
        "reference": "See POWER_MARKERS",
        "creatine_response": {
            "CC": "RR - Good creatine responder",
            "TT": "XX - May be non-responder to creatine"
        },
        "actionable": {
            "TT": [
                "May be creatine non-responder (XX genotype)",
                "Higher endogenous creatine stores possible",
                "Can still try creatine - individual variation exists"
            ],
            "CC": [
                "Likely creatine responder",
                "Standard loading/maintenance protocols effective"
            ]
        }
    },
    "rs7118114": {
        "gene": "CNDP1",
        "variant": "Carnosinase D18S880",
        "function": "Carnosine degradation - affects muscle buffering",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35},
        "effect": "Associated with creatine and beta-alanine response",
        "category": "supplement_response",
        "evidence": "moderate",
        "pmid": ["20303006"]
    },
}

# =============================================================================
# CAFFEINE ERGOGENIC
# =============================================================================

CAFFEINE_ERGOGENIC_MARKERS = {
    "rs762551": {
        "gene": "CYP1A2",
        "variant": "*1F",
        "function": "Caffeine metabolism",
        "risk_allele": "C",
        "frequency": {"EUR": 0.32},
        "effect": {
            "AA": "Fast metabolizer - CAFFEINE IMPROVES performance",
            "AC": "Intermediate - Modest benefit",
            "CC": "Slow metabolizer - Caffeine may IMPAIR performance"
        },
        "category": "caffeine_ergogenic",
        "evidence": "strong",
        "pmid": ["16522833", "26219105", "28770561"],
        "actionable": {
            "AA": [
                "CAFFEINE ERGOGENIC - improves exercise performance",
                "3-6 mg/kg caffeine 30-60min pre-exercise",
                "Can use for endurance and power activities",
                "Timing less critical (fast clearance)"
            ],
            "CC": [
                "Caffeine may NOT improve or may IMPAIR performance",
                "Slow metabolism means longer CNS stimulation",
                "Consider avoiding pre-competition or using very low dose",
                "May interfere with recovery if used post-exercise"
            ]
        }
    },
    "rs5751876": {
        "gene": "ADORA2A",
        "variant": "1976T>C",
        "function": "Adenosine receptor - caffeine target",
        "risk_allele": "T",
        "frequency": {"EUR": 0.40},
        "effect": {
            "TT": "High caffeine sensitivity - may get anxiety from pre-workout caffeine",
            "CC": "Normal sensitivity - standard ergogenic effect"
        },
        "category": "caffeine_ergogenic",
        "evidence": "moderate",
        "pmid": ["18088379", "22038822"],
        "actionable": {
            "TT": [
                "High caffeine sensitivity",
                "May experience anxiety/jitters that impair performance",
                "Lower doses or avoid altogether for competition"
            ]
        }
    },
}

# =============================================================================
# HEAT TOLERANCE
# =============================================================================

HEAT_MARKERS = {
    "rs1043618": {
        "gene": "HSPA1A",
        "variant": "HSP70-1 (G190C)",
        "function": "Heat shock protein 70 - thermal stress response",
        "risk_allele": "C",
        "frequency": {"EUR": 0.35, "AFR": 0.20},
        "effect": {
            "CC": "Reduced heat shock protein response",
            "GG": "Better thermal stress adaptation"
        },
        "category": "heat",
        "evidence": "moderate",
        "pmid": ["21418717", "24627588"],
        "actionable": {
            "CC": [
                "May have reduced heat tolerance",
                "Extra heat acclimatization important",
                "Careful hydration in hot conditions",
                "Consider pre-cooling strategies"
            ]
        }
    },
    "rs6457452": {
        "gene": "HSPA1B",
        "variant": "HSP70-2",
        "function": "Heat shock protein response",
        "risk_allele": "C",
        "frequency": {"EUR": 0.25},
        "effect": "Associated with heat tolerance",
        "category": "heat",
        "evidence": "moderate",
        "pmid": ["24627588"]
    },
}

# =============================================================================
# ALTITUDE ADAPTATION
# =============================================================================

ALTITUDE_MARKERS = {
    "rs1867785": {
        "gene": "EPAS1",
        "variant": "HIF-2α (Endothelial PAS domain protein 1)",
        "function": "Hypoxia response - EPO regulation",
        "risk_allele": "G",
        "frequency": {"EUR": 0.20, "TIB": 0.80},  # Very high in Tibetans
        "effect": {
            "GG": "Tibetan-like adaptation - better high-altitude tolerance"
        },
        "category": "altitude",
        "evidence": "definitive",
        "pmid": ["20466884", "24974814"],
        "actionable": {
            "GG": [
                "May have better altitude adaptation",
                "Could respond well to altitude training",
                "Less altitude sickness risk"
            ]
        }
    },
    "rs479200": {
        "gene": "EGLN1",
        "variant": "PHD2 (Prolyl hydroxylase)",
        "function": "HIF degradation - oxygen sensing",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15, "AND": 0.50},  # Elevated in Andean populations
        "effect": "Associated with altitude adaptation",
        "category": "altitude",
        "evidence": "strong",
        "pmid": ["20466884", "17997608"]
    },
    "rs11549465": {
        "gene": "HIF1A",
        "reference": "See ENDURANCE_MARKERS",
        "altitude_effect": "Pro582Ser associated with altitude performance"
    },
}

# =============================================================================
# CONCUSSION / BRAIN INJURY
# =============================================================================

CONCUSSION_MARKERS = {
    "rs429358": {
        "gene": "APOE",
        "variant": "APOE ε4 marker (C112R)",
        "function": "Apolipoprotein E - brain repair and recovery",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15, "AFR": 0.25, "EAS": 0.08},
        "effect": {
            "ε4 carriers": "Slower recovery from TBI, higher long-term risk"
        },
        "category": "concussion",
        "evidence": "strong",
        "pmid": ["9069288", "25818573", "27445615"],
        "actionable": {
            "e4_carrier": [
                "SLOWER concussion recovery",
                "Extended rest periods may be needed after head injury",
                "Consider avoiding high-concussion-risk sports",
                "More conservative return-to-play protocols",
                "Higher CTE risk with repeated head impacts"
            ]
        }
    },
    "rs1800497": {
        "gene": "ANKK1/DRD2",
        "variant": "Taq1A",
        "function": "Dopamine receptor - affects recovery",
        "risk_allele": "A",
        "frequency": {"EUR": 0.20, "EAS": 0.40},
        "effect": {
            "AA": "A1/A1 - May have prolonged symptoms post-concussion"
        },
        "category": "concussion",
        "evidence": "moderate",
        "pmid": ["25818573"],
        "actionable": {
            "AA": [
                "May have prolonged post-concussion symptoms",
                "Cognitive symptoms may persist longer"
            ]
        }
    },
}

# =============================================================================
# HYDRATION
# =============================================================================

HYDRATION_MARKERS = {
    "rs4340": {
        "gene": "ACE",
        "reference": "See ENDURANCE_MARKERS",
        "hydration_effect": {
            "DD": "May need more sodium - salt-sensitive",
            "II": "More efficient fluid regulation"
        }
    },
    "rs3741559": {
        "gene": "AQP1",
        "variant": "Aquaporin 1",
        "function": "Water channel protein",
        "risk_allele": "A",
        "frequency": {"EUR": 0.15},
        "effect": "Associated with exercise-induced hyponatremia risk",
        "category": "hydration",
        "evidence": "moderate",
        "pmid": ["25688896"],
        "actionable": {
            "AA": [
                "May be at higher risk of exercise-induced hyponatremia",
                "Don't overhydrate during long events",
                "Use electrolytes, not just water"
            ]
        }
    },
}

# =============================================================================
# COMBINE ALL ATHLETIC MARKERS
# =============================================================================

ATHLETIC_COMPREHENSIVE_MARKERS = {
    **POWER_MARKERS,
    **ENDURANCE_MARKERS,
    **LACTATE_MARKERS,
    **RECOVERY_MARKERS,
    **TENDON_MARKERS,
    **BONE_MARKERS,
    **CREATINE_MARKERS,
    **CAFFEINE_ERGOGENIC_MARKERS,
    **HEAT_MARKERS,
    **ALTITUDE_MARKERS,
    **CONCUSSION_MARKERS,
    **HYDRATION_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def calculate_athletic_profile(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate comprehensive athletic genetic profile."""
    
    power_score = 0
    endurance_score = 0
    
    # ACTN3 - major effect
    actn3 = genotypes.get("rs1815739", "CT")
    if actn3 == "CC":  # RR
        power_score += 3
    elif actn3 == "TT":  # XX
        endurance_score += 3
    else:
        power_score += 1
        endurance_score += 1
    
    # ACE I/D
    ace = genotypes.get("rs4340")
    if ace == "II":
        endurance_score += 2
    elif ace == "DD":
        power_score += 2
    else:
        power_score += 1
        endurance_score += 1
    
    # PPARGC1A
    pgc1a = genotypes.get("rs8192678", "GA")
    if pgc1a == "GG":
        endurance_score += 1
    
    # Determine profile
    total = power_score + endurance_score
    ratio = power_score / max(total, 1)
    
    if ratio >= 0.7:
        profile = AthleticProfile.POWER_DOMINANT
    elif ratio >= 0.55:
        profile = AthleticProfile.POWER_ENDURANCE
    elif ratio >= 0.45:
        profile = AthleticProfile.BALANCED
    elif ratio >= 0.3:
        profile = AthleticProfile.ENDURANCE_POWER
    else:
        profile = AthleticProfile.ENDURANCE_DOMINANT
    
    # Sport recommendations
    sport_map = {
        AthleticProfile.POWER_DOMINANT: ["Sprint events", "Weightlifting", "Shot put", "Football (positions requiring burst)"],
        AthleticProfile.POWER_ENDURANCE: ["Swimming", "Soccer", "Basketball", "Martial arts", "400m"],
        AthleticProfile.BALANCED: ["Tennis", "Triathlon", "CrossFit", "Soccer", "Most sports"],
        AthleticProfile.ENDURANCE_POWER: ["Middle distance", "Cycling", "Rowing", "Soccer midfielder"],
        AthleticProfile.ENDURANCE_DOMINANT: ["Marathon", "Ultrarunning", "Cycling TT", "Triathlon", "XC skiing"],
    }
    
    return {
        "power_score": power_score,
        "endurance_score": endurance_score,
        "profile": profile.value,
        "recommended_sports": sport_map.get(profile, []),
        "actn3_status": "RR (power)" if actn3 == "CC" else "XX (endurance)" if actn3 == "TT" else "RX (mixed)",
        "ace_status": "II (endurance)" if ace == "II" else "DD (power)" if ace == "DD" else "ID (mixed)" if ace else "Unknown"
    }

def calculate_injury_risk(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate soft tissue injury risk."""
    
    risk_score = 0
    risk_factors = []
    
    # COL1A1
    col1a1 = genotypes.get("rs1800012", "GG")
    if col1a1 == "TT":
        risk_score += 3
        risk_factors.append("COL1A1 TT - High collagen injury risk")
    elif col1a1 == "GT":
        risk_score += 1
        risk_factors.append("COL1A1 GT - Moderate risk")
    
    # COL5A1
    col5a1 = genotypes.get("rs12722", "CC")
    if col5a1 == "TT":
        risk_score += 2
        risk_factors.append("COL5A1 TT - ACL/Achilles risk elevated")
    elif col5a1 == "CT":
        risk_score += 1
    
    # TNC
    tnc = genotypes.get("rs13321", "CC")
    if tnc == "TT":
        risk_score += 1.5
        risk_factors.append("TNC TT - Achilles tendinopathy risk")
    
    # Determine level
    if risk_score >= 5:
        level = InjuryRiskLevel.HIGH
    elif risk_score >= 3:
        level = InjuryRiskLevel.ELEVATED
    elif risk_score >= 1.5:
        level = InjuryRiskLevel.MODERATE
    else:
        level = InjuryRiskLevel.LOW
    
    return {
        "risk_score": round(risk_score, 1),
        "risk_level": level.value,
        "risk_factors": risk_factors,
        "recommendations": [
            "Extended warm-up protocols" if risk_score >= 3 else None,
            "Collagen supplementation (10-15g + vitamin C)" if risk_score >= 3 else None,
            "Progressive loading essential" if risk_score >= 3 else None,
            "Consider prehabilitation exercises" if risk_score >= 1.5 else None,
            "ACL prevention program recommended" if col5a1 in ["TT", "CT"] else None,
            "Eccentric calf exercises for Achilles protection" if tnc == "TT" else None,
        ],
        "pmid": ["16477526", "20460586"]
    }

def calculate_recovery_profile(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate recovery speed profile."""
    
    il6 = genotypes.get("rs1800795", "GC")
    tnf = genotypes.get("rs1800629", "GG")
    
    if il6 == "CC" or tnf in ["AA", "GA"]:
        speed = RecoverySpeed.SLOW
        rest_days = "48-72+ hours between intense sessions"
    elif il6 == "GG" and tnf == "GG":
        speed = RecoverySpeed.FAST
        rest_days = "24-48 hours, can train more frequently"
    else:
        speed = RecoverySpeed.NORMAL
        rest_days = "24-48 hours standard recovery"
    
    return {
        "recovery_speed": speed.value,
        "recommended_rest": rest_days,
        "inflammation_profile": "High" if speed == RecoverySpeed.SLOW else "Normal" if speed == RecoverySpeed.NORMAL else "Low",
        "recommendations": [
            "Prioritize sleep (8+ hours)" if speed == RecoverySpeed.SLOW else "Standard sleep",
            "Anti-inflammatory nutrition (omega-3, tart cherry)" if speed == RecoverySpeed.SLOW else None,
            "Cold water immersion may help" if speed == RecoverySpeed.SLOW else None,
            "Can tolerate higher training frequency" if speed == RecoverySpeed.FAST else None,
        ],
        "pmid": ["15616363", "23632419"]
    }

def generate_athletic_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate comprehensive athletic genetics report."""
    
    profile = calculate_athletic_profile(genotypes)
    injury_risk = calculate_injury_risk(genotypes)
    recovery = calculate_recovery_profile(genotypes)
    
    # Caffeine ergogenic
    cyp1a2 = genotypes.get("rs762551", "AC")
    caffeine_helps = cyp1a2 == "AA"
    
    # Creatine response
    actn3 = genotypes.get("rs1815739", "CT")
    creatine_responder = actn3 != "TT"
    
    return {
        "athletic_profile": profile,
        "injury_risk": injury_risk,
        "recovery_profile": recovery,
        "supplement_responses": {
            "caffeine_ergogenic": caffeine_helps,
            "caffeine_note": "Improves performance" if caffeine_helps else "May not help or could impair",
            "creatine_responder": creatine_responder,
            "creatine_note": "Likely responder" if creatine_responder else "May be non-responder (XX genotype)"
        },
        "markers_analyzed": sum(1 for rs in ATHLETIC_COMPREHENSIVE_MARKERS if rs in genotypes),
    }

# Export
__all__ = [
    'ATHLETIC_COMPREHENSIVE_MARKERS',
    'POWER_MARKERS',
    'ENDURANCE_MARKERS',
    'LACTATE_MARKERS',
    'RECOVERY_MARKERS',
    'TENDON_MARKERS',
    'BONE_MARKERS',
    'CREATINE_MARKERS',
    'CAFFEINE_ERGOGENIC_MARKERS',
    'HEAT_MARKERS',
    'ALTITUDE_MARKERS',
    'CONCUSSION_MARKERS',
    'HYDRATION_MARKERS',
    'AthleticProfile',
    'RecoverySpeed',
    'InjuryRiskLevel',
    'calculate_athletic_profile',
    'calculate_injury_risk',
    'calculate_recovery_profile',
    'generate_athletic_report',
]
