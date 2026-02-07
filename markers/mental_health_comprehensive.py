"""
Comprehensive Mental Health Genetics v5.0

Complete coverage of:
- Depression (SLC6A4/5-HTTLPR, BDNF, FKBP5, CRHR1)
- Antidepressant response (CYP2D6, CYP2C19, SLC6A4, HTR2A)
- Anxiety (COMT, NPSR1, GABRA2, RGS2)
- ADHD (DRD4, DAT1/SLC6A3, SNAP25, COMT)
- Bipolar (CACNA1C, ANK3, BDNF)
- Stress response (COMT, FKBP5, NR3C1, CRHR1)
- Resilience (NPY, OXTR, COMT)
- Social bonding/empathy (OXTR, CD38, AVPR1A)
- Aggression/impulsivity (MAOA, 5-HTTLPR, COMT)
- Personality traits (DRD4, 5-HTTLPR, ANKK1)
- Addiction susceptibility

All markers with PMID references and clinical considerations.

IMPORTANT DISCLAIMER: Psychiatric genetics are complex. These markers 
contribute small effects. Mental health is influenced by environment, 
experiences, and many genes. Use for awareness, not diagnosis.
"""

from enum import Enum
from typing import Dict, List, Any, Optional

class RiskLevel(Enum):
    ELEVATED = "elevated"
    MODERATE = "moderate"
    AVERAGE = "average"
    REDUCED = "reduced"

class StressResponseType(Enum):
    WARRIOR = "warrior"      # Val/Val COMT - stress resilient, less focused
    WORRIER = "worrier"      # Met/Met COMT - stress sensitive, detail-oriented
    INTERMEDIATE = "intermediate"

class EvidenceStrength(Enum):
    REPLICATED = "replicated"      # Multiple large studies
    MODERATE = "moderate"          # Some replication
    PRELIMINARY = "preliminary"    # Single/small studies
    CONTROVERSIAL = "controversial" # Mixed results

# =============================================================================
# DISCLAIMER
# =============================================================================
MENTAL_HEALTH_DISCLAIMER = """
IMPORTANT: Psychiatric genetics are highly complex and polygenic.
- Individual variants have SMALL effects (typically <1% variance explained)
- Environment, life experiences, and epigenetics play major roles
- These results are NOT diagnostic
- Do not use to predict or label mental health conditions
- Discuss with healthcare provider if concerned
- Useful for: medication selection, self-understanding, awareness
"""

# =============================================================================
# DEPRESSION
# =============================================================================

DEPRESSION_MARKERS = {
    "rs25531": {
        "gene": "SLC6A4",
        "variant": "5-HTTLPR L/S equivalent (A/G SNP)",
        "function": "Serotonin transporter expression",
        "risk_allele": "A",  # Approximates S allele
        "frequency": {"EUR": 0.45, "AFR": 0.25, "EAS": 0.75},
        "effect": {
            "AA": "S/S equivalent - lower serotonin transporter, higher depression risk with stress",
            "AG": "L/S equivalent - intermediate",
            "GG": "L/L equivalent - higher serotonin transporter"
        },
        "category": "depression",
        "evidence": EvidenceStrength.CONTROVERSIAL,
        "pmid": ["12869766", "14993431", "28855158"],
        "note": "Gene x environment interaction - depression risk WITH stressful life events",
        "actionable": {
            "AA": [
                "May be more sensitive to stressful life events",
                "Stress management especially important",
                "May respond well to SSRIs (more transporter to block)",
                "Mindfulness and CBT particularly beneficial",
                "Build strong social support networks"
            ]
        }
    },
    "rs6265": {
        "gene": "BDNF",
        "variant": "Val66Met",
        "function": "Brain-derived neurotrophic factor - neuroplasticity",
        "risk_allele": "T",  # Met allele
        "frequency": {"EUR": 0.20, "EAS": 0.45, "AFR": 0.02},
        "effect": {
            "TT": "Met/Met - Reduced BDNF secretion, poorer stress resilience",
            "CT": "Val/Met - Intermediate",
            "CC": "Val/Val - Normal BDNF function"
        },
        "category": "depression",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["12802784", "15341763", "16689191"],
        "actionable": {
            "TT": [
                "Lower activity-dependent BDNF release",
                "Exercise is especially important (increases BDNF)",
                "May have poorer episodic memory",
                "Antidepressant response may differ",
                "Chronic stress particularly harmful"
            ],
            "CT": [
                "Intermediate BDNF function",
                "Exercise still beneficial for BDNF"
            ]
        }
    },
    "rs1360780": {
        "gene": "FKBP5",
        "variant": "FKBP5 stress response",
        "function": "Glucocorticoid receptor sensitivity - stress axis",
        "risk_allele": "T",
        "frequency": {"EUR": 0.30, "AFR": 0.10},
        "effect": {
            "TT": "Prolonged cortisol response, higher PTSD/depression risk with trauma",
            "CC": "Normal stress axis regulation"
        },
        "category": "depression",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["15289279", "21362375", "23831920"],
        "note": "Gene x environment: risk with childhood adversity",
        "actionable": {
            "TT": [
                "Prolonged cortisol response to stress",
                "Higher risk of PTSD and depression with childhood trauma",
                "Stress management and therapy especially important",
                "Early intervention after trauma beneficial"
            ]
        }
    },
    "rs242924": {
        "gene": "CRHR1",
        "variant": "Corticotropin-releasing hormone receptor 1",
        "function": "HPA axis regulation",
        "risk_allele": "T",
        "frequency": {"EUR": 0.45},
        "effect": "Protective haplotype against depression with childhood maltreatment",
        "category": "depression",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["16650898", "18701695"]
    },
}

# =============================================================================
# ANTIDEPRESSANT RESPONSE
# =============================================================================

ANTIDEPRESSANT_MARKERS = {
    "rs4244285": {
        "gene": "CYP2C19",
        "variant": "*2",
        "function": "SSRI metabolism",
        "risk_allele": "A",
        "frequency": {"EUR": 0.15, "EAS": 0.30},
        "effect": {
            "AA": "Poor metabolizer - higher SSRI levels, more side effects",
            "AG": "Intermediate metabolizer",
            "GG": "Normal metabolizer"
        },
        "category": "antidepressant",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["25974703", "23486447"],
        "drugs_affected": ["citalopram", "escitalopram", "sertraline", "amitriptyline"],
        "actionable": {
            "AA": [
                "Poor CYP2C19 metabolizer",
                "Lower doses of SSRIs needed (especially escitalopram)",
                "Higher risk of side effects at standard doses",
                "Citalopram max 20mg/day for CYP2C19 PMs"
            ]
        }
    },
    "rs3892097": {
        "gene": "CYP2D6",
        "variant": "*4",
        "function": "Tricyclic and some SSRI metabolism",
        "risk_allele": "A",
        "frequency": {"EUR": 0.20},
        "effect": {
            "AA": "Poor metabolizer - reduced metabolism of TCAs, some SSRIs"
        },
        "category": "antidepressant",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["27997040", "28520587"],
        "drugs_affected": ["nortriptyline", "desipramine", "venlafaxine", "fluoxetine", "paroxetine"],
        "actionable": {
            "AA": [
                "CYP2D6 poor metabolizer",
                "Tricyclic antidepressants: use 50% initial dose",
                "Venlafaxine: may have reduced conversion to active metabolite"
            ]
        }
    },
    "rs6313": {
        "gene": "HTR2A",
        "variant": "T102C (5-HT2A receptor)",
        "function": "Serotonin 2A receptor density",
        "risk_allele": "T",
        "frequency": {"EUR": 0.40, "EAS": 0.60},
        "effect": {
            "TT": "Higher 5-HT2A receptor expression - may have better SSRI response",
            "CC": "Lower receptor density"
        },
        "category": "antidepressant",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["12140781", "16855965"],
        "actionable": {
            "TT": [
                "May have better response to SSRIs",
                "Higher serotonin receptor density"
            ]
        }
    },
    "rs25531": {
        "gene": "SLC6A4",
        "reference": "See DEPRESSION_MARKERS",
        "antidepressant_effect": {
            "AA": "May respond well to SSRIs (more target)",
            "GG": "May need higher SSRI doses or different class"
        }
    },
}

# =============================================================================
# ANXIETY
# =============================================================================

ANXIETY_MARKERS = {
    "rs4680": {
        "gene": "COMT",
        "variant": "Val158Met",
        "function": "Catecholamine degradation in prefrontal cortex",
        "risk_allele": "A",  # Met = slow = worrier
        "frequency": {"EUR": 0.48, "EAS": 0.30, "AFR": 0.40},
        "effect": {
            "AA": "Met/Met - WORRIER: Higher PFC dopamine, better focus but anxiety-prone",
            "GG": "Val/Val - WARRIOR: Lower PFC dopamine, stress resilient but less focused",
            "AG": "Intermediate phenotype"
        },
        "category": "anxiety",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["12716966", "17008817", "14517761"],
        "stress_type": {
            "AA": StressResponseType.WORRIER,
            "AG": StressResponseType.INTERMEDIATE,
            "GG": StressResponseType.WARRIOR
        },
        "actionable": {
            "AA": [
                "WORRIER phenotype - high prefrontal dopamine",
                "Better working memory and attention",
                "BUT more anxiety-prone under stress",
                "Stress management essential",
                "May be MORE sensitive to stimulants/caffeine",
                "Mindfulness particularly helpful",
                "May need LESS dopamine-boosting interventions"
            ],
            "GG": [
                "WARRIOR phenotype - stress resilient",
                "May need more stimulation to focus",
                "Better performance under pressure",
                "May benefit from dopamine-supporting supplements",
                "Lower anxiety baseline"
            ]
        }
    },
    "rs324981": {
        "gene": "NPSR1",
        "variant": "Asn107Ile",
        "function": "Neuropeptide S receptor - arousal/anxiety",
        "risk_allele": "T",  # Ile allele
        "frequency": {"EUR": 0.45},
        "effect": {
            "TT": "Ile/Ile - Higher anxiety trait, panic disorder risk",
            "CC": "Asn/Asn - Lower anxiety"
        },
        "category": "anxiety",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["19088741", "16648366"],
        "actionable": {
            "TT": [
                "May have higher anxiety trait",
                "Elevated panic disorder risk",
                "Anxiety management techniques important"
            ]
        }
    },
    "rs279858": {
        "gene": "GABRA2",
        "variant": "GABA-A receptor alpha-2 subunit",
        "function": "GABAergic inhibition",
        "risk_allele": "G",
        "frequency": {"EUR": 0.40, "AFR": 0.20},
        "effect": {
            "GG": "Associated with anxiety and alcohol dependence"
        },
        "category": "anxiety",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["15051850", "17005050"],
        "note": "Also associated with alcoholism risk"
    },
    "rs4606": {
        "gene": "RGS2",
        "variant": "Regulator of G-protein signaling 2",
        "function": "Anxiety signaling modulation",
        "risk_allele": "C",
        "frequency": {"EUR": 0.50},
        "effect": {
            "CC": "Associated with higher anxiety trait"
        },
        "category": "anxiety",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["18358332", "19996040"]
    },
}

# =============================================================================
# ADHD
# =============================================================================

ADHD_MARKERS = {
    "rs1800955": {
        "gene": "DRD4",
        "variant": "DRD4 -521 C/T",
        "function": "Dopamine D4 receptor promoter",
        "risk_allele": "T",
        "frequency": {"EUR": 0.45, "EAS": 0.30},
        "effect": {
            "TT": "Lower DRD4 expression - associated with novelty seeking and ADHD"
        },
        "category": "adhd",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["15858146", "10945462"],
        "note": "7-repeat VNTR (not available on chip) has stronger ADHD association"
    },
    "rs28363170": {
        "gene": "SLC6A3",
        "variant": "DAT1 3' UTR VNTR proxy",
        "function": "Dopamine transporter expression",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30},
        "effect": "Associated with ADHD and stimulant response",
        "category": "adhd",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["11923833", "15534621"],
        "note": "10-repeat VNTR associated with ADHD"
    },
    "rs3746544": {
        "gene": "SNAP25",
        "variant": "SNAP-25 (synaptosomal protein)",
        "function": "Neurotransmitter release",
        "risk_allele": "T",
        "frequency": {"EUR": 0.40},
        "effect": {
            "TT": "Associated with ADHD, affects synaptic transmission"
        },
        "category": "adhd",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["12823878", "19557953"]
    },
    "rs4680": {
        "gene": "COMT",
        "reference": "See ANXIETY_MARKERS",
        "adhd_effect": {
            "AA": "Met/Met - Better attention but anxiety risk",
            "GG": "Val/Val - May benefit more from stimulant medications"
        },
        "actionable": {
            "adhd_context": [
                "Val/Val may respond better to stimulant ADHD medications",
                "Met/Met may be more sensitive to stimulant side effects"
            ]
        }
    },
}

# =============================================================================
# BIPOLAR DISORDER
# =============================================================================

BIPOLAR_MARKERS = {
    "rs1006737": {
        "gene": "CACNA1C",
        "variant": "L-type calcium channel alpha-1C",
        "function": "Neuronal calcium signaling",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35, "EAS": 0.10},
        "effect": {
            "AA": "Associated with bipolar disorder and schizophrenia",
            "note": "One of strongest GWAS hits for bipolar"
        },
        "category": "bipolar",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["18349090", "21926974"],
        "actionable": {
            "AA": [
                "Elevated bipolar disorder genetic risk marker",
                "NOT diagnostic - many with this genotype don't develop bipolar",
                "Be aware of mood cycling patterns",
                "Calcium channel blockers being studied as treatment"
            ]
        }
    },
    "rs10994336": {
        "gene": "ANK3",
        "variant": "Ankyrin G",
        "function": "Neuronal structure and ion channel clustering",
        "risk_allele": "T",
        "frequency": {"EUR": 0.10},
        "effect": "Strong GWAS hit for bipolar disorder",
        "category": "bipolar",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["18711365", "21926974"]
    },
    "rs6265": {
        "gene": "BDNF",
        "reference": "See DEPRESSION_MARKERS",
        "bipolar_effect": "Met allele associated with earlier onset bipolar"
    },
}

# =============================================================================
# STRESS RESPONSE / RESILIENCE
# =============================================================================

STRESS_MARKERS = {
    "rs4680": {
        "gene": "COMT",
        "reference": "See ANXIETY_MARKERS (primary reference)"
    },
    "rs1360780": {
        "gene": "FKBP5",
        "reference": "See DEPRESSION_MARKERS (primary reference)"
    },
    "rs41423247": {
        "gene": "NR3C1",
        "variant": "Glucocorticoid receptor BclI",
        "function": "Cortisol receptor sensitivity",
        "risk_allele": "G",
        "frequency": {"EUR": 0.35},
        "effect": {
            "GG": "Increased GR sensitivity - metabolic effects, stress response"
        },
        "category": "stress",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["16950492", "17360155"]
    },
    "rs16147": {
        "gene": "NPY",
        "variant": "Neuropeptide Y promoter",
        "function": "Stress resilience, anxiety modulation",
        "risk_allele": "T",
        "frequency": {"EUR": 0.45},
        "effect": {
            "TT": "Lower NPY expression - reduced stress resilience",
            "CC": "Higher NPY - protective against stress"
        },
        "category": "stress",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["21795697", "21106944"],
        "actionable": {
            "CC": [
                "Higher neuropeptide Y levels",
                "Associated with stress resilience",
                "Protective factor"
            ]
        }
    },
}

# =============================================================================
# SOCIAL / EMPATHY (OXYTOCIN SYSTEM)
# =============================================================================

SOCIAL_MARKERS = {
    "rs53576": {
        "gene": "OXTR",
        "variant": "Oxytocin receptor A/G",
        "function": "Oxytocin receptor expression and function",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35, "EAS": 0.60},
        "effect": {
            "GG": "Higher empathy, better social cognition, BUT also more stress response",
            "AA": "Lower empathy trait, less socially sensitive"
        },
        "category": "social",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["19015103", "21595907", "21908899"],
        "actionable": {
            "GG": [
                "Higher trait empathy and emotional sensitivity",
                "Better social cognition and reading emotions",
                "May also be more affected by others' stress",
                "Strong social support more important"
            ],
            "AA": [
                "Lower baseline empathy trait",
                "Less affected by social stress",
                "May need to actively work on emotional attunement"
            ]
        }
    },
    "rs3796863": {
        "gene": "CD38",
        "variant": "CD38 oxytocin release",
        "function": "Regulates oxytocin release",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30},
        "effect": {
            "AA": "Lower oxytocin release, associated with autism risk"
        },
        "category": "social",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["20824728", "22422981"]
    },
    "rs7632287": {
        "gene": "AVPR1A",
        "variant": "Vasopressin receptor 1A",
        "function": "Social bonding, pair bonding",
        "risk_allele": "G",
        "frequency": {"EUR": 0.20},
        "effect": "Associated with pair bonding and relationship satisfaction",
        "category": "social",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["18776892", "20566715"]
    },
}

# =============================================================================
# AGGRESSION / IMPULSIVITY
# =============================================================================

AGGRESSION_MARKERS = {
    "rs909525": {
        "gene": "MAOA",
        "variant": "MAOA upstream SNP",
        "function": "Monoamine oxidase A - serotonin/dopamine breakdown",
        "risk_allele": "T",
        "frequency": {"EUR": 0.40},
        "effect": {
            "note": "MAOA VNTR (not on chip) is the key variant for 'warrior gene'",
            "low_activity": "Low MAOA activity + childhood maltreatment → aggression risk"
        },
        "category": "aggression",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["12161658", "17189274"],
        "note": "Gene x environment - only with childhood maltreatment"
    },
    "rs25531": {
        "gene": "SLC6A4",
        "reference": "See DEPRESSION_MARKERS",
        "aggression_effect": "S allele + stress → impulsivity"
    },
    "rs4680": {
        "gene": "COMT",
        "reference": "See ANXIETY_MARKERS",
        "aggression_effect": {
            "GG": "Val/Val may have higher aggression under alcohol"
        }
    },
}

# =============================================================================
# PERSONALITY TRAITS
# =============================================================================

PERSONALITY_MARKERS = {
    "rs1800955": {
        "gene": "DRD4",
        "reference": "See ADHD_MARKERS",
        "personality_effect": {
            "novelty_seeking": "T allele associated with novelty seeking",
            "note": "7-repeat VNTR is the classic novelty-seeking variant"
        }
    },
    "rs25531": {
        "gene": "SLC6A4",
        "reference": "See DEPRESSION_MARKERS",
        "personality_effect": {
            "AA": "S/S - Higher neuroticism, harm avoidance"
        }
    },
    "rs1800497": {
        "gene": "ANKK1/DRD2",
        "variant": "Taq1A",
        "function": "Dopamine D2 receptor density",
        "risk_allele": "A",  # A1 allele
        "frequency": {"EUR": 0.20, "EAS": 0.40, "AFR": 0.35},
        "effect": {
            "AA": "A1/A1 - Lower D2 receptor density, reward deficiency",
            "AG": "A1/A2 - Intermediate",
            "GG": "A2/A2 - Normal D2 density"
        },
        "category": "personality",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["8807664", "2906084"],
        "actionable": {
            "AA": [
                "Lower dopamine receptor density",
                "May seek more stimulation/rewards",
                "Higher addiction susceptibility",
                "May need more intense stimuli to feel rewarded"
            ]
        }
    },
}

# =============================================================================
# ADDICTION SUSCEPTIBILITY
# =============================================================================

ADDICTION_MARKERS = {
    # Alcohol
    "rs671": {
        "gene": "ALDH2",
        "reference": "See nutrition_comprehensive.py",
        "addiction_effect": {
            "GA/AA": "PROTECTIVE against alcoholism (flushing)",
            "GG": "No protection"
        }
    },
    "rs279858": {
        "gene": "GABRA2",
        "reference": "See ANXIETY_MARKERS",
        "addiction_effect": {
            "GG": "Associated with alcohol dependence"
        }
    },
    "rs1799971": {
        "gene": "OPRM1",
        "variant": "A118G (μ-opioid receptor)",
        "function": "Opioid and reward response",
        "risk_allele": "G",
        "frequency": {"EUR": 0.15, "EAS": 0.40},
        "effect": {
            "GG": "Higher alcohol reward, may drink more for effect",
            "AG": "Better naltrexone response for alcohol use disorder"
        },
        "category": "addiction",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["21412232", "17329694"],
        "actionable": {
            "AG/GG": [
                "G allele: May have stronger alcohol reward",
                "BUT also better response to naltrexone for AUD",
                "If struggling with alcohol, naltrexone may help"
            ]
        }
    },
    # Nicotine
    "rs16969968": {
        "gene": "CHRNA5",
        "variant": "D398N (nicotinic receptor)",
        "function": "Nicotinic acetylcholine receptor",
        "risk_allele": "A",
        "frequency": {"EUR": 0.35, "EAS": 0.05, "AFR": 0.05},
        "effect": {
            "AA": "Higher nicotine dependence risk, heavier smoking"
        },
        "category": "addiction",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["18385739", "21559498"],
        "actionable": {
            "AA": [
                "Higher nicotine addiction susceptibility",
                "If you smoke - higher dependence risk",
                "May need more intensive cessation support",
                "Varenicline (Chantix) may be particularly effective"
            ]
        }
    },
    "rs4105144": {
        "gene": "CYP2A6",
        "variant": "Nicotine metabolism",
        "function": "Nicotine breakdown speed",
        "risk_allele": "C",
        "frequency": {"EUR": 0.05, "EAS": 0.15},
        "effect": {
            "CC": "Slow nicotine metabolism - may smoke less",
            "TT": "Fast metabolism - more cigarettes needed"
        },
        "category": "addiction",
        "evidence": EvidenceStrength.REPLICATED,
        "pmid": ["16845390", "21559498"]
    },
    # Opioid
    "rs1799971_opioid": {
        "gene": "OPRM1",
        "reference": "See above - same variant",
        "opioid_effect": {
            "GG": "Reduced opioid analgesia, may need higher doses",
            "AA": "Normal opioid response"
        }
    },
    # Cannabis
    "rs1049353": {
        "gene": "CNR1",
        "variant": "Cannabinoid receptor 1",
        "function": "Endocannabinoid signaling",
        "risk_allele": "A",
        "frequency": {"EUR": 0.25},
        "effect": "Associated with cannabis dependence susceptibility",
        "category": "addiction",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["12072389", "14557543"]
    },
    "rs324420": {
        "gene": "FAAH",
        "variant": "Fatty acid amide hydrolase C385A",
        "function": "Endocannabinoid breakdown",
        "risk_allele": "A",  # Reduced FAAH
        "frequency": {"EUR": 0.20, "AFR": 0.45},
        "effect": {
            "AA": "Lower FAAH - higher endocannabinoid levels, lower anxiety",
            "note": "May be less drawn to cannabis (already have high endocannabinoids)"
        },
        "category": "addiction",
        "evidence": EvidenceStrength.MODERATE,
        "pmid": ["20966023", "19474099"],
        "actionable": {
            "AA": [
                "Higher baseline endocannabinoid levels",
                "May have lower anxiety naturally",
                "Paradoxically may be less drawn to cannabis"
            ]
        }
    },
    # General addiction
    "rs1800497": {
        "gene": "ANKK1/DRD2",
        "reference": "See PERSONALITY_MARKERS",
        "addiction_effect": {
            "AA": "A1/A1 - General addiction vulnerability (reward deficiency)"
        }
    },
}

# =============================================================================
# COMBINE ALL MENTAL HEALTH MARKERS
# =============================================================================

MENTAL_HEALTH_COMPREHENSIVE_MARKERS = {
    **DEPRESSION_MARKERS,
    **ANTIDEPRESSANT_MARKERS,
    **ANXIETY_MARKERS,
    **ADHD_MARKERS,
    **BIPOLAR_MARKERS,
    **STRESS_MARKERS,
    **SOCIAL_MARKERS,
    **AGGRESSION_MARKERS,
    **PERSONALITY_MARKERS,
    **ADDICTION_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def determine_stress_type(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Determine warrior vs worrier phenotype from COMT."""
    comt = genotypes.get("rs4680", "AG")
    
    if comt == "AA":
        stress_type = StressResponseType.WORRIER
        description = "Worrier - Higher prefrontal dopamine, detail-oriented, anxiety-prone"
    elif comt == "GG":
        stress_type = StressResponseType.WARRIOR
        description = "Warrior - Lower prefrontal dopamine, stress resilient, may need more stimulation"
    else:
        stress_type = StressResponseType.INTERMEDIATE
        description = "Intermediate - Balanced stress response"
    
    return {
        "stress_type": stress_type.value,
        "comt_genotype": comt,
        "description": description,
        "recommendations": [
            "Stress management important" if stress_type == StressResponseType.WORRIER else "May thrive under pressure",
            "May be caffeine sensitive" if stress_type == StressResponseType.WORRIER else "May tolerate stimulants well",
            "Detail-oriented work suits well" if stress_type == StressResponseType.WORRIER else "High-pressure environments may suit"
        ],
        "pmid": ["12716966", "17008817"]
    }

def assess_depression_vulnerability(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Assess genetic depression vulnerability markers."""
    
    factors = []
    risk_score = 0
    
    # 5-HTTLPR
    sert = genotypes.get("rs25531", "AG")
    if sert == "AA":
        factors.append("5-HTTLPR S/S equivalent - stress sensitivity")
        risk_score += 1
    
    # BDNF
    bdnf = genotypes.get("rs6265", "CC")
    if bdnf == "TT":
        factors.append("BDNF Met/Met - reduced neuroplasticity")
        risk_score += 1
    elif bdnf == "CT":
        factors.append("BDNF Val/Met - mild effect")
        risk_score += 0.5
    
    # FKBP5
    fkbp5 = genotypes.get("rs1360780", "CC")
    if fkbp5 == "TT":
        factors.append("FKBP5 risk variant - prolonged stress response")
        risk_score += 1
    
    # Determine level
    if risk_score >= 2:
        level = RiskLevel.ELEVATED
    elif risk_score >= 1:
        level = RiskLevel.MODERATE
    else:
        level = RiskLevel.AVERAGE
    
    return {
        "risk_level": level.value,
        "risk_factors": factors,
        "protective_recommendations": [
            "Strong social support is protective",
            "Regular exercise boosts BDNF",
            "Stress management and mindfulness",
            "Early intervention after trauma",
            "Therapy (CBT/DBT) particularly effective"
        ],
        "disclaimer": MENTAL_HEALTH_DISCLAIMER
    }

def assess_addiction_vulnerability(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Assess addiction susceptibility markers."""
    
    factors = []
    
    # DRD2
    drd2 = genotypes.get("rs1800497", "GG")
    if drd2 in ["AA", "AG"]:
        factors.append("DRD2 A1 allele - reward deficiency, general addiction risk")
    
    # GABRA2
    gabra2 = genotypes.get("rs279858", "AA")
    if gabra2 == "GG":
        factors.append("GABRA2 - alcohol dependence susceptibility")
    
    # CHRNA5
    chrna5 = genotypes.get("rs16969968", "GG")
    if chrna5 == "AA":
        factors.append("CHRNA5 - nicotine dependence susceptibility")
    
    # OPRM1
    oprm1 = genotypes.get("rs1799971", "AA")
    if oprm1 in ["GG", "AG"]:
        factors.append("OPRM1 G allele - may have stronger alcohol reward")
    
    return {
        "vulnerability_factors": factors,
        "substances_of_concern": [
            "Alcohol" if gabra2 == "GG" or oprm1 in ["GG", "AG"] else None,
            "Nicotine" if chrna5 == "AA" else None,
            "General reward-seeking" if drd2 in ["AA", "AG"] else None,
        ],
        "protective_strategies": [
            "Awareness of vulnerability",
            "Avoid or minimize high-risk substances",
            "Seek help early if developing dependence",
            "Naltrexone may help for alcohol (OPRM1 G carriers)"
        ],
        "disclaimer": MENTAL_HEALTH_DISCLAIMER
    }

def generate_mental_health_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate comprehensive mental health genetics report."""
    
    stress_type = determine_stress_type(genotypes)
    depression_risk = assess_depression_vulnerability(genotypes)
    addiction_risk = assess_addiction_vulnerability(genotypes)
    
    # Social/empathy
    oxtr = genotypes.get("rs53576", "AG")
    empathy_level = "high" if oxtr == "GG" else "low" if oxtr == "AA" else "moderate"
    
    return {
        "disclaimer": MENTAL_HEALTH_DISCLAIMER,
        "stress_phenotype": stress_type,
        "depression_vulnerability": depression_risk,
        "addiction_vulnerability": addiction_risk,
        "social_profile": {
            "empathy_trait": empathy_level,
            "oxtr_genotype": oxtr
        },
        "markers_analyzed": sum(1 for rs in MENTAL_HEALTH_COMPREHENSIVE_MARKERS if rs in genotypes),
        "important_note": "These are risk factors, not diagnoses. Environment and experiences matter enormously."
    }

# Export
__all__ = [
    'MENTAL_HEALTH_COMPREHENSIVE_MARKERS',
    'MENTAL_HEALTH_DISCLAIMER',
    'DEPRESSION_MARKERS',
    'ANTIDEPRESSANT_MARKERS',
    'ANXIETY_MARKERS',
    'ADHD_MARKERS',
    'BIPOLAR_MARKERS',
    'STRESS_MARKERS',
    'SOCIAL_MARKERS',
    'AGGRESSION_MARKERS',
    'PERSONALITY_MARKERS',
    'ADDICTION_MARKERS',
    'RiskLevel',
    'StressResponseType',
    'EvidenceStrength',
    'determine_stress_type',
    'assess_depression_vulnerability',
    'assess_addiction_vulnerability',
    'generate_mental_health_report',
]
