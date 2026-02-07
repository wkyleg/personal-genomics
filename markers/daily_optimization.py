"""
Daily Optimization Genetics v5.0

Comprehensive coverage of:
- Chronotype genetics (CLOCK, PER1, PER2, PER3, CRY1, CRY2, ARNTL)
- Complete caffeine metabolism (CYP1A2, ADORA2A, AHR)
- Exercise timing optimization
- Meal timing / circadian metabolism

Generates personalized daily schedule recommendations.
"""

from enum import Enum
from typing import Dict, List, Any, Optional
from dataclasses import dataclass

class Chronotype(Enum):
    DEFINITE_MORNING = "definite_morning"    # 5:30-6:30am wake
    MODERATE_MORNING = "moderate_morning"    # 6:30-7:30am wake
    INTERMEDIATE = "intermediate"            # 7:30-8:30am wake
    MODERATE_EVENING = "moderate_evening"    # 8:30-10am wake
    DEFINITE_EVENING = "definite_evening"    # 10am+ natural wake

class CaffeineMetabolism(Enum):
    ULTRAFAST = "ultrafast"      # <2 hours half-life
    FAST = "fast"                # 2-4 hours half-life
    NORMAL = "normal"            # 4-6 hours half-life
    SLOW = "slow"                # 6-9 hours half-life
    ULTRASLOW = "ultraslow"      # >9 hours half-life

class ExerciseTimingType(Enum):
    MORNING_OPTIMAL = "morning_optimal"
    AFTERNOON_OPTIMAL = "afternoon_optimal"
    EVENING_OPTIMAL = "evening_optimal"
    FLEXIBLE = "flexible"

# =============================================================================
# CHRONOTYPE / CIRCADIAN GENES
# =============================================================================

CHRONOTYPE_MARKERS = {
    # CLOCK - Master circadian clock gene
    "rs1801260": {
        "gene": "CLOCK",
        "variant": "3111T>C",
        "function": "Master circadian transcription factor",
        "risk_allele": "C",
        "frequency": {"EUR": 0.28, "EAS": 0.15, "AFR": 0.20},
        "effect": {
            "CC": "Evening preference, delayed sleep phase",
            "TC": "Moderate evening tendency",
            "TT": "Morning preference"
        },
        "category": "chronotype",
        "evidence": "strong",
        "pmid": ["11927384", "17550340", "23657362"],
        "chronotype_effect": {
            "CC": -2,  # Strong evening
            "TC": -1,  # Moderate evening
            "TT": 0    # Neutral/morning
        },
        "actionable": {
            "CC": [
                "Natural evening type (night owl)",
                "May function best with later schedule",
                "Important decisions/creative work better in evening",
                "Exercise in late afternoon/evening optimal",
                "Light exposure in morning helps shift earlier"
            ],
            "TT": [
                "Natural morning type (early bird)",
                "Peak alertness in morning hours",
                "Schedule important work for morning"
            ]
        }
    },
    "rs12649507": {
        "gene": "CLOCK",
        "variant": "CLOCK GWAS hit",
        "function": "Circadian rhythm regulation",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30},
        "effect": "Morning preference marker",
        "category": "chronotype",
        "evidence": "definitive",
        "pmid": ["30804565", "26835600"],
        "chronotype_effect": {
            "AA": 1,
            "AG": 0.5,
            "GG": 0
        }
    },
    
    # PER1 - Period Circadian Protein 1
    "rs2735611": {
        "gene": "PER1",
        "variant": "PER1 polymorphism",
        "function": "Circadian rhythm period length",
        "risk_allele": "G",
        "frequency": {"EUR": 0.40, "EAS": 0.30},
        "effect": "Affects morning/evening preference",
        "category": "chronotype",
        "evidence": "moderate",
        "pmid": ["19487505", "15927956"],
        "chronotype_effect": {
            "GG": 0.5,
            "AG": 0,
            "AA": -0.5
        }
    },
    
    # PER2 - Period Circadian Protein 2
    "rs35333999": {
        "gene": "PER2",
        "variant": "PER2 variant",
        "function": "Core clock component - period determination",
        "risk_allele": "C",
        "frequency": {"EUR": 0.20},
        "effect": {
            "CC": "Associated with advanced sleep phase",
            "note": "PER2 mutations cause familial advanced sleep phase disorder"
        },
        "category": "chronotype",
        "evidence": "strong",
        "pmid": ["11232563", "17535983"],
        "chronotype_effect": {
            "CC": 1.5,  # Strong morning
            "CT": 0.5,
            "TT": 0
        },
        "actionable": {
            "CC": [
                "Strong morning preference (natural early riser)",
                "May wake naturally at 4-5am",
                "Peak performance in early morning",
                "Evening light exposure can help delay"
            ]
        }
    },
    
    # PER3 - Period Circadian Protein 3
    "rs228697": {
        "gene": "PER3",
        "variant": "VNTR-associated SNP",
        "function": "Sleep homeostasis and circadian timing",
        "risk_allele": "C",
        "frequency": {"EUR": 0.08, "AFR": 0.02},
        "effect": {
            "CC": "Short PER3 VNTR linked - evening type, needs less sleep",
            "GG": "Long VNTR linked - morning type, more sleep pressure"
        },
        "category": "chronotype",
        "evidence": "strong",
        "pmid": ["12841365", "16596785"],
        "chronotype_effect": {
            "CC": -1,
            "CG": -0.5,
            "GG": 1
        },
        "actionable": {
            "CC": [
                "May function well on less sleep",
                "Evening type tendency",
                "Less impacted by sleep deprivation"
            ],
            "GG": [
                "Higher sleep need (8+ hours)",
                "Morning type",
                "More vulnerable to sleep deprivation effects"
            ]
        }
    },
    "rs228729": {
        "gene": "PER3",
        "variant": "PER3 expression variant",
        "function": "Circadian rhythm amplitude",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35},
        "effect": "Affects chronotype",
        "category": "chronotype",
        "evidence": "moderate",
        "pmid": ["30804565"],
        "chronotype_effect": {
            "AA": 0.5,
            "AG": 0,
            "GG": -0.5
        }
    },
    
    # CRY1 - Cryptochrome 1
    "rs2287161": {
        "gene": "CRY1",
        "variant": "CRY1 variant",
        "function": "Negative feedback loop of circadian clock",
        "risk_allele": "C",
        "frequency": {"EUR": 0.35, "EAS": 0.25},
        "effect": {
            "CC": "Delayed sleep phase tendency",
            "note": "CRY1 mutation causes DSPD"
        },
        "category": "chronotype",
        "evidence": "strong",
        "pmid": ["28351990", "30804565"],
        "chronotype_effect": {
            "CC": -1.5,
            "CG": -0.5,
            "GG": 0
        },
        "actionable": {
            "CC": [
                "Strong evening preference / delayed sleep phase",
                "May have difficulty falling asleep at conventional times",
                "Light therapy in morning can help",
                "Melatonin 2-3 hours before desired sleep may help"
            ]
        }
    },
    
    # CRY2 - Cryptochrome 2
    "rs10838524": {
        "gene": "CRY2",
        "variant": "CRY2 variant",
        "function": "Circadian rhythm regulation",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30},
        "effect": "Affects chronotype and sleep timing",
        "category": "chronotype",
        "evidence": "moderate",
        "pmid": ["30804565", "27046545"],
        "chronotype_effect": {
            "AA": -0.5,
            "AG": 0,
            "GG": 0.5
        }
    },
    
    # ARNTL (BMAL1) - Core clock transcription factor
    "rs7950226": {
        "gene": "ARNTL",
        "variant": "BMAL1 variant",
        "function": "Positive arm of circadian clock",
        "risk_allele": "A",
        "frequency": {"EUR": 0.40},
        "effect": "Affects chronotype and metabolic rhythm",
        "category": "chronotype",
        "evidence": "moderate",
        "pmid": ["30804565", "16547163"],
        "chronotype_effect": {
            "AA": 0.5,
            "AG": 0,
            "GG": -0.5
        }
    },
    "rs4757144": {
        "gene": "ARNTL",
        "variant": "BMAL1 metabolic variant",
        "function": "Circadian metabolism regulation",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35},
        "effect": "Affects metabolic timing",
        "category": "chronotype",
        "evidence": "moderate",
        "pmid": ["16547163"],
        "metabolic_effect": "May affect glucose tolerance timing"
    },
    
    # Additional GWAS chronotype loci
    "rs73598374": {
        "gene": "ADA (adenosine deaminase)",
        "variant": "Asp8Asn",
        "function": "Adenosine metabolism - sleep pressure",
        "risk_allele": "T",
        "frequency": {"EUR": 0.05},
        "effect": {
            "TT": "Reduced adenosine breakdown, higher sleep pressure",
            "note": "Rare homozygotes have extreme morning preference"
        },
        "category": "chronotype",
        "evidence": "strong",
        "pmid": ["30804565", "15964771"],
        "chronotype_effect": {
            "TT": 2,  # Very strong morning
            "CT": 0.5,
            "CC": 0
        },
        "actionable": {
            "TT": [
                "Very strong morning type",
                "High sleep pressure (deep sleeper)",
                "May need more sleep than average"
            ]
        }
    },
}

# =============================================================================
# CAFFEINE METABOLISM
# =============================================================================

CAFFEINE_MARKERS = {
    # CYP1A2 - Primary caffeine metabolizer
    "rs762551": {
        "gene": "CYP1A2",
        "variant": "*1F (intron 1)",
        "function": "Liver caffeine metabolism (95% of caffeine)",
        "risk_allele": "C",
        "frequency": {"EUR": 0.32, "EAS": 0.25, "AFR": 0.40},
        "effect": {
            "AA": "Fast metabolizer (induced by smoking, charred foods)",
            "AC": "Intermediate metabolizer",
            "CC": "Slow metabolizer"
        },
        "category": "caffeine",
        "evidence": "definitive",
        "pmid": ["17132150", "16522833", "20485253"],
        "caffeine_effect": {
            "AA": CaffeineMetabolism.FAST,
            "AC": CaffeineMetabolism.NORMAL,
            "CC": CaffeineMetabolism.SLOW
        },
        "actionable": {
            "AA": [
                "Fast caffeine metabolizer (half-life ~3-4 hours)",
                "Can tolerate caffeine later in day without sleep disruption",
                "May need more frequent caffeine for sustained effect",
                "Moderate coffee consumption may have cardiovascular benefit",
                "Note: smoking induces CYP1A2 making it even faster"
            ],
            "CC": [
                "Slow caffeine metabolizer (half-life ~6-9 hours)",
                "CAFFEINE CUTOFF: No later than 12-2pm for good sleep",
                "Higher anxiety/jitteriness with caffeine",
                "Increased cardiovascular risk with >2 cups/day",
                "Coffee afternoon = likely sleep disruption"
            ],
            "AC": [
                "Intermediate caffeine metabolizer",
                "Caffeine cutoff: 2-4pm",
                "Moderate consumption generally fine"
            ]
        }
    },
    "rs2069514": {
        "gene": "CYP1A2",
        "variant": "*1C",
        "function": "CYP1A2 expression",
        "risk_allele": "A",
        "frequency": {"EAS": 0.25, "SAS": 0.10},
        "effect": {
            "AA": "Decreased CYP1A2 expression - slower metabolism"
        },
        "category": "caffeine",
        "evidence": "moderate",
        "pmid": ["8476986", "17132150"],
        "caffeine_effect": {
            "AA": CaffeineMetabolism.SLOW
        }
    },
    "rs12720461": {
        "gene": "CYP1A2",
        "variant": "*1K",
        "function": "CYP1A2 expression",
        "risk_allele": "G",
        "frequency": {"EUR": 0.02, "AFR": 0.03},
        "effect": "Decreased expression - slower metabolism",
        "category": "caffeine",
        "evidence": "moderate",
        "pmid": ["17132150"]
    },
    
    # ADORA2A - Adenosine A2A receptor
    "rs5751876": {
        "gene": "ADORA2A",
        "variant": "1976T>C (c.1083)",
        "function": "Adenosine receptor - caffeine binding site",
        "risk_allele": "T",
        "frequency": {"EUR": 0.40, "EAS": 0.30, "AFR": 0.25},
        "effect": {
            "TT": "Higher caffeine sensitivity (anxiety, insomnia)",
            "CT": "Moderate sensitivity",
            "CC": "Normal caffeine response"
        },
        "category": "caffeine",
        "evidence": "strong",
        "pmid": ["18088379", "22038822", "21857820"],
        "sensitivity_effect": {
            "TT": "high_sensitivity",
            "CT": "moderate_sensitivity",
            "CC": "normal_sensitivity"
        },
        "actionable": {
            "TT": [
                "HIGH caffeine sensitivity - affects receptor binding, not metabolism",
                "More anxiety/nervousness with caffeine",
                "Greater sleep disruption even with morning caffeine",
                "May need to limit total caffeine (<100mg/day)",
                "Consider decaf or tea instead of coffee"
            ],
            "CC": [
                "Normal caffeine sensitivity",
                "Can typically enjoy caffeine without excessive anxiety"
            ]
        }
    },
    "rs2298383": {
        "gene": "ADORA2A",
        "variant": "ADORA2A variant",
        "function": "Adenosine receptor expression",
        "risk_allele": "C",
        "frequency": {"EUR": 0.40},
        "effect": "Modulates caffeine response",
        "category": "caffeine",
        "evidence": "moderate",
        "pmid": ["22038822"]
    },
    "rs3761422": {
        "gene": "ADORA2A",
        "variant": "Anxiety-associated variant",
        "function": "Adenosine receptor",
        "risk_allele": "T",
        "frequency": {"EUR": 0.35},
        "effect": "Associated with caffeine-induced anxiety",
        "category": "caffeine",
        "evidence": "moderate",
        "pmid": ["21857820"]
    },
    
    # AHR - Aryl hydrocarbon receptor
    "rs2066853": {
        "gene": "AHR",
        "variant": "R554K",
        "function": "Transcription factor regulating CYP1A2",
        "risk_allele": "A",
        "frequency": {"EUR": 0.10, "EAS": 0.30},
        "effect": {
            "AA": "Reduced CYP1A2 induction - slower caffeine clearance"
        },
        "category": "caffeine",
        "evidence": "moderate",
        "pmid": ["19103647", "20932375"],
        "caffeine_effect": {
            "AA": CaffeineMetabolism.SLOW
        },
        "actionable": {
            "AA": [
                "Reduced ability to induce CYP1A2",
                "May have slower caffeine metabolism",
                "Diet/smoking has less effect on metabolism"
            ]
        }
    },
}

# =============================================================================
# EXERCISE TIMING
# =============================================================================

EXERCISE_TIMING_MARKERS = {
    "rs1801260": {
        "gene": "CLOCK",
        "reference": "See CHRONOTYPE_MARKERS",
        "exercise_timing": {
            "CC": "Evening exercise optimal (5-8pm)",
            "TT": "Morning exercise optimal (6-10am)",
            "TC": "Flexible, afternoon may be best"
        }
    },
    "rs12649507": {
        "gene": "CLOCK",
        "reference": "See CHRONOTYPE_MARKERS",
        "exercise_timing": {
            "AA": "Morning exercise optimal",
            "GG": "Afternoon/evening exercise may be better"
        }
    },
    
    # Performance circadian markers
    "rs1144566": {
        "gene": "PER2",
        "variant": "PER2 performance variant",
        "function": "Circadian rhythm and muscle performance",
        "risk_allele": "G",
        "frequency": {"EUR": 0.35},
        "effect": "Affects timing of peak physical performance",
        "category": "exercise_timing",
        "evidence": "moderate",
        "pmid": ["17535983", "26835600"],
        "exercise_timing": {
            "GG": "Morning peak performance",
            "AA": "Evening peak performance",
            "AG": "Afternoon peak"
        }
    },
    
    # Muscle clock genes
    "rs4757144": {
        "gene": "ARNTL",
        "variant": "BMAL1 muscle variant",
        "function": "Muscle circadian clock",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35},
        "effect": "Affects muscle performance timing",
        "category": "exercise_timing",
        "evidence": "moderate",
        "pmid": ["16547163"]
    },
}

# =============================================================================
# MEAL TIMING / CIRCADIAN METABOLISM
# =============================================================================

MEAL_TIMING_MARKERS = {
    "rs4757144": {
        "gene": "ARNTL",
        "variant": "BMAL1 metabolic",
        "function": "Metabolic circadian regulation",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35},
        "effect": {
            "AA": "May have better glucose tolerance in morning",
            "GG": "More flexible glucose handling timing"
        },
        "category": "meal_timing",
        "evidence": "moderate",
        "pmid": ["16547163", "22675159"],
        "actionable": {
            "AA": [
                "Glucose tolerance better in morning",
                "Front-load calories earlier in day",
                "Avoid large late-night meals"
            ]
        }
    },
    "rs1801260": {
        "gene": "CLOCK",
        "variant": "CLOCK 3111T>C",
        "function": "Metabolic timing",
        "risk_allele": "C",
        "frequency": {"EUR": 0.28},
        "effect": {
            "CC": "Associated with later eating times, weight gain",
            "note": "Evening types with CC more prone to metabolic issues"
        },
        "category": "meal_timing",
        "evidence": "strong",
        "pmid": ["22675159", "23277449"],
        "actionable": {
            "CC": [
                "Natural late eating tendency - associated with weight gain",
                "Try to eat dinner by 7pm despite evening preference",
                "Eating window earlier in day may help metabolic health",
                "More vulnerable to negative effects of shift work"
            ]
        }
    },
    "rs2943641": {
        "gene": "IRS1",
        "variant": "Insulin receptor substrate",
        "function": "Insulin signaling and timing",
        "risk_allele": "C",
        "frequency": {"EUR": 0.60},
        "effect": {
            "CC": "Better insulin sensitivity in morning",
            "TT": "More consistent insulin response across day"
        },
        "category": "meal_timing",
        "evidence": "moderate",
        "pmid": ["22231480"],
        "actionable": {
            "CC": [
                "Insulin sensitivity peaks in morning",
                "Carbohydrates better tolerated at breakfast/lunch",
                "Lighter, lower-carb dinners recommended"
            ]
        }
    },
    # Melatonin receptor - affects glucose
    "rs10830963": {
        "gene": "MTNR1B",
        "variant": "Melatonin receptor 1B",
        "function": "Melatonin signaling in pancreas",
        "risk_allele": "G",
        "frequency": {"EUR": 0.30, "EAS": 0.45, "AFR": 0.10},
        "effect": {
            "GG": "Impaired insulin secretion, especially at night",
            "note": "MTNR1B G allele = T2D risk + late-eating problems"
        },
        "category": "meal_timing",
        "evidence": "definitive",
        "pmid": ["19060907", "25681352", "26391390"],
        "actionable": {
            "GG": [
                "AVOID late-night eating - especially carbs",
                "Melatonin impairs your glucose tolerance more than others",
                "Finish eating 3-4 hours before bed",
                "Morning/early afternoon carbs better",
                "Higher T2D risk if late eater"
            ],
            "AG": [
                "Moderate caution with late-night carbs",
                "Earlier eating window beneficial"
            ]
        }
    },
}

# =============================================================================
# COMBINE ALL DAILY OPTIMIZATION MARKERS
# =============================================================================

DAILY_OPTIMIZATION_MARKERS = {
    **CHRONOTYPE_MARKERS,
    **CAFFEINE_MARKERS,
    **EXERCISE_TIMING_MARKERS,
    **MEAL_TIMING_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def calculate_chronotype_score(genotypes: Dict[str, str]) -> tuple[float, Chronotype]:
    """
    Calculate chronotype score from genetic markers.
    Negative = evening type, Positive = morning type.
    Returns (score, chronotype_enum).
    """
    score = 0.0
    markers_found = 0
    
    chronotype_effects = {
        "rs1801260": {"CC": -2, "TC": -1, "TT": 0},
        "rs12649507": {"AA": 1, "AG": 0.5, "GG": 0},
        "rs2735611": {"GG": 0.5, "AG": 0, "AA": -0.5},
        "rs35333999": {"CC": 1.5, "CT": 0.5, "TT": 0},
        "rs228697": {"CC": -1, "CG": -0.5, "GG": 1},
        "rs228729": {"AA": 0.5, "AG": 0, "GG": -0.5},
        "rs2287161": {"CC": -1.5, "CG": -0.5, "GG": 0},
        "rs10838524": {"AA": -0.5, "AG": 0, "GG": 0.5},
        "rs7950226": {"AA": 0.5, "AG": 0, "GG": -0.5},
        "rs73598374": {"TT": 2, "CT": 0.5, "CC": 0},
    }
    
    for rsid, effects in chronotype_effects.items():
        geno = genotypes.get(rsid)
        if geno and geno in effects:
            score += effects[geno]
            markers_found += 1
    
    # Normalize if few markers
    if markers_found > 0:
        score = score * (10 / max(markers_found, 5))  # Scale to expected range
    
    # Determine chronotype
    if score >= 3:
        chronotype = Chronotype.DEFINITE_MORNING
    elif score >= 1:
        chronotype = Chronotype.MODERATE_MORNING
    elif score >= -1:
        chronotype = Chronotype.INTERMEDIATE
    elif score >= -3:
        chronotype = Chronotype.MODERATE_EVENING
    else:
        chronotype = Chronotype.DEFINITE_EVENING
    
    return round(score, 1), chronotype

def calculate_caffeine_profile(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate complete caffeine response profile."""
    
    # Metabolism (CYP1A2)
    cyp1a2_main = genotypes.get("rs762551", "AC")
    cyp1a2_1c = genotypes.get("rs2069514", "GG")
    
    if cyp1a2_main == "CC" or cyp1a2_1c == "AA":
        metabolism = CaffeineMetabolism.SLOW
        half_life_hours = "6-9"
        cutoff_time = "12:00 PM (noon)"
    elif cyp1a2_main == "AA" and cyp1a2_1c != "AA":
        metabolism = CaffeineMetabolism.FAST
        half_life_hours = "3-4"
        cutoff_time = "4:00-6:00 PM"
    else:
        metabolism = CaffeineMetabolism.NORMAL
        half_life_hours = "4-6"
        cutoff_time = "2:00-4:00 PM"
    
    # Sensitivity (ADORA2A)
    adora2a = genotypes.get("rs5751876", "CC")
    
    if adora2a == "TT":
        sensitivity = "high"
        anxiety_risk = "elevated"
        max_daily = "100-200mg (1-2 small cups)"
    elif adora2a == "CT":
        sensitivity = "moderate"
        anxiety_risk = "moderate"
        max_daily = "200-400mg (2-4 cups)"
    else:
        sensitivity = "normal"
        anxiety_risk = "normal"
        max_daily = "400mg (4 cups) - standard limit"
    
    # CVD risk consideration
    cvd_risk = "elevated" if metabolism == CaffeineMetabolism.SLOW and sensitivity == "normal" else "normal"
    
    return {
        "metabolism": metabolism.value,
        "half_life_hours": half_life_hours,
        "sensitivity": sensitivity,
        "anxiety_risk": anxiety_risk,
        "caffeine_cutoff_time": cutoff_time,
        "recommended_max_daily": max_daily,
        "cardiovascular_risk": cvd_risk,
        "recommendations": [
            f"Caffeine metabolism: {metabolism.value} (half-life ~{half_life_hours} hours)",
            f"Sensitivity: {sensitivity}",
            f"Last caffeine by: {cutoff_time}",
            f"Recommended max: {max_daily}",
            "Slow metabolizers: >2 cups/day may increase cardiovascular risk" if metabolism == CaffeineMetabolism.SLOW else "",
            "High sensitivity: May experience anxiety/jitters even at low doses" if sensitivity == "high" else "",
        ],
        "pmid": ["17132150", "16522833", "18088379"]
    }

def calculate_optimal_exercise_time(genotypes: Dict[str, str], chronotype: Chronotype) -> Dict[str, Any]:
    """Determine optimal exercise timing based on genetics."""
    
    clock_3111 = genotypes.get("rs1801260", "TC")
    
    # Map chronotype to exercise timing
    if chronotype in [Chronotype.DEFINITE_MORNING, Chronotype.MODERATE_MORNING]:
        optimal_window = "6:00 AM - 10:00 AM"
        peak_performance = "8:00 AM - 9:00 AM"
        timing_type = ExerciseTimingType.MORNING_OPTIMAL
    elif chronotype in [Chronotype.DEFINITE_EVENING, Chronotype.MODERATE_EVENING]:
        optimal_window = "4:00 PM - 8:00 PM"
        peak_performance = "5:00 PM - 7:00 PM"
        timing_type = ExerciseTimingType.EVENING_OPTIMAL
    else:
        optimal_window = "2:00 PM - 6:00 PM"
        peak_performance = "3:00 PM - 5:00 PM"
        timing_type = ExerciseTimingType.AFTERNOON_OPTIMAL
    
    # Adjust for CLOCK variant
    if clock_3111 == "CC":
        # Evening CLOCK - shift recommendations later
        notes = "CLOCK CC genotype: Your body clock runs late. Evening exercise aligns with your circadian peak."
    elif clock_3111 == "TT":
        notes = "CLOCK TT genotype: Your body clock favors morning. Early exercise aligns with your rhythm."
    else:
        notes = "Intermediate chronotype: Flexible exercise timing, but avoid late evening for sleep quality."
    
    return {
        "optimal_window": optimal_window,
        "peak_performance_time": peak_performance,
        "timing_type": timing_type.value,
        "notes": notes,
        "recommendations": [
            "Body temperature peaks ~5-7pm for most people (strength/power)",
            "Morning exercise: Better adherence, cortisol utilization",
            "Evening exercise: Better performance, but may disrupt sleep if too late",
            f"For your chronotype ({chronotype.value}): {optimal_window} is genetically optimal"
        ]
    }

def calculate_optimal_meal_timing(genotypes: Dict[str, str], chronotype: Chronotype) -> Dict[str, Any]:
    """Determine optimal meal timing based on genetics."""
    
    mtnr1b = genotypes.get("rs10830963", "CC")
    clock = genotypes.get("rs1801260", "TC")
    irs1 = genotypes.get("rs2943641", "CT")
    
    warnings = []
    
    # MTNR1B is critical
    if mtnr1b == "GG":
        late_eating_risk = "HIGH"
        dinner_cutoff = "6:00 PM"
        carb_timing = "Morning/early afternoon ONLY"
        warnings.append("⚠️ MTNR1B GG: Late eating significantly impairs glucose tolerance")
    elif mtnr1b == "AG":
        late_eating_risk = "Moderate"
        dinner_cutoff = "7:00 PM"
        carb_timing = "Earlier in day preferred"
    else:
        late_eating_risk = "Normal"
        dinner_cutoff = "8:00 PM"
        carb_timing = "Flexible"
    
    # IRS1 consideration
    if irs1 == "CC":
        carb_timing = "Morning/lunch best" if carb_timing == "Flexible" else carb_timing
    
    # CLOCK evening type
    if clock == "CC":
        warnings.append("CLOCK CC: Natural late eater - discipline required for metabolic health")
    
    return {
        "late_eating_risk": late_eating_risk,
        "dinner_cutoff": dinner_cutoff,
        "carbohydrate_timing": carb_timing,
        "recommended_eating_window": "7:00 AM - 7:00 PM" if mtnr1b == "GG" else "7:00 AM - 8:00 PM",
        "warnings": warnings,
        "recommendations": [
            "Largest meal at breakfast or lunch (not dinner)",
            f"Finish eating by {dinner_cutoff}",
            f"Carbohydrate timing: {carb_timing}",
            "Protein can be more evenly distributed",
            "Fast 12+ hours overnight for metabolic health"
        ],
        "pmid": ["25681352", "26391390", "22675159"]
    }

def generate_daily_optimization_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate comprehensive daily optimization report."""
    
    chrono_score, chronotype = calculate_chronotype_score(genotypes)
    caffeine_profile = calculate_caffeine_profile(genotypes)
    exercise_timing = calculate_optimal_exercise_time(genotypes, chronotype)
    meal_timing = calculate_optimal_meal_timing(genotypes, chronotype)
    
    # Generate ideal schedule
    if chronotype in [Chronotype.DEFINITE_MORNING, Chronotype.MODERATE_MORNING]:
        ideal_wake = "5:30 - 6:30 AM"
        ideal_sleep = "9:30 - 10:30 PM"
        peak_cognitive = "8:00 AM - 12:00 PM"
    elif chronotype in [Chronotype.DEFINITE_EVENING, Chronotype.MODERATE_EVENING]:
        ideal_wake = "8:30 - 10:00 AM"
        ideal_sleep = "12:00 - 1:30 AM"
        peak_cognitive = "4:00 PM - 10:00 PM"
    else:
        ideal_wake = "7:00 - 8:00 AM"
        ideal_sleep = "11:00 PM - 12:00 AM"
        peak_cognitive = "10:00 AM - 2:00 PM"
    
    return {
        "chronotype": {
            "genetic_score": chrono_score,
            "type": chronotype.value,
            "interpretation": f"Genetic {chronotype.value.replace('_', ' ')} chronotype",
            "ideal_wake_time": ideal_wake,
            "ideal_sleep_time": ideal_sleep,
            "peak_cognitive_window": peak_cognitive,
        },
        "caffeine": caffeine_profile,
        "exercise_timing": exercise_timing,
        "meal_timing": meal_timing,
        "daily_schedule": {
            "wake": ideal_wake,
            "peak_cognitive": peak_cognitive,
            "exercise": exercise_timing["optimal_window"],
            "dinner_by": meal_timing["dinner_cutoff"],
            "caffeine_cutoff": caffeine_profile["caffeine_cutoff_time"],
            "sleep": ideal_sleep,
        },
        "markers_analyzed": sum(1 for rs in DAILY_OPTIMIZATION_MARKERS if rs in genotypes),
    }

# Export
__all__ = [
    'DAILY_OPTIMIZATION_MARKERS',
    'CHRONOTYPE_MARKERS',
    'CAFFEINE_MARKERS',
    'EXERCISE_TIMING_MARKERS',
    'MEAL_TIMING_MARKERS',
    'Chronotype',
    'CaffeineMetabolism',
    'ExerciseTimingType',
    'calculate_chronotype_score',
    'calculate_caffeine_profile',
    'calculate_optimal_exercise_time',
    'calculate_optimal_meal_timing',
    'generate_daily_optimization_report',
]
