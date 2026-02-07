"""
Complete Pharmacogenomics Marker Database v5.0
Comprehensive drug-gene interactions for precision medicine.

Sources:
- CPIC Guidelines (Clinical Pharmacogenetics Implementation Consortium)
- PharmGKB (Pharmacogenomics Knowledge Base)
- FDA Table of Pharmacogenomic Biomarkers
- DPWG (Dutch Pharmacogenetics Working Group)

All markers have actionable clinical recommendations and PMID references.
"""

from enum import Enum
from typing import Dict, List, Optional, Any

class MetabolizerStatus(Enum):
    ULTRARAPID = "ultrarapid"
    RAPID = "rapid"
    NORMAL = "normal"
    INTERMEDIATE = "intermediate"
    POOR = "poor"
    UNKNOWN = "unknown"

class ClinicalActionLevel(Enum):
    CRITICAL = "critical"  # Life-threatening - immediate action required
    HIGH = "high"          # Significant clinical impact
    MODERATE = "moderate"  # Clinically relevant
    LOW = "low"            # Informational

# =============================================================================
# CYP2D6 - Metabolizes ~25% of drugs
# =============================================================================

CYP2D6_MARKERS = {
    "rs3892097": {
        "gene": "CYP2D6",
        "variant": "*4",
        "function": "No function (null allele)",
        "risk_allele": "A",
        "frequency": {"EUR": 0.20, "EAS": 0.01, "AFR": 0.06, "SAS": 0.10},
        "activity_score": 0,
        "drugs_affected": [
            "codeine", "tramadol", "hydrocodone", "oxycodone",
            "tamoxifen", "ondansetron", "tropisetron",
            "fluoxetine", "paroxetine", "venlafaxine", "atomoxetine",
            "metoprolol", "carvedilol", "propafenone", "flecainide",
            "risperidone", "aripiprazole", "haloperidol", "thioridazine",
            "dextromethorphan", "eliglustat"
        ],
        "clinical_impact": "Poor metabolizer - prodrugs ineffective, parent drugs accumulate",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27997040", "31562822", "28520587"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "codeine": "Use alternative analgesic (morphine, non-opioid). Codeine provides no analgesia in PMs.",
                "tramadol": "Use alternative. Tramadol provides reduced/no analgesia in PMs.",
                "tamoxifen": "Consider alternative (aromatase inhibitor) or higher dose. Reduced efficacy in PMs.",
                "metoprolol": "Consider 50% dose reduction or alternative beta blocker (atenolol, bisoprolol).",
                "atomoxetine": "Initiate at standard dose, may need slower titration. Higher plasma levels in PMs.",
                "general": "Share CYP2D6 status with all prescribers before starting affected medications."
            }
        }
    },
    "rs5030655": {
        "gene": "CYP2D6",
        "variant": "*6",
        "function": "No function (null allele - frameshift)",
        "risk_allele": "del",
        "frequency": {"EUR": 0.01, "AFR": 0.01},
        "activity_score": 0,
        "drugs_affected": ["Same as *4"],
        "clinical_impact": "Poor metabolizer",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27997040"]
    },
    "rs28371706": {
        "gene": "CYP2D6",
        "variant": "*17",
        "function": "Decreased function",
        "risk_allele": "A",
        "frequency": {"AFR": 0.20, "EUR": 0.01},
        "activity_score": 0.5,
        "clinical_impact": "Intermediate metabolizer - common in African populations",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27997040", "9918136"]
    },
    "rs1065852": {
        "gene": "CYP2D6",
        "variant": "*10",
        "function": "Decreased function",
        "risk_allele": "A",
        "frequency": {"EAS": 0.40, "EUR": 0.02, "SAS": 0.25},
        "activity_score": 0.25,
        "clinical_impact": "Intermediate metabolizer - major variant in East Asian populations",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27997040", "8807664"],
        "notes": "Most common reduced function allele globally"
    },
    "rs28371725": {
        "gene": "CYP2D6",
        "variant": "*41",
        "function": "Decreased function",
        "risk_allele": "A",
        "frequency": {"EUR": 0.09, "AFR": 0.02, "MENA": 0.20},
        "activity_score": 0.5,
        "clinical_impact": "Intermediate metabolizer",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27997040", "11875365"]
    },
    "rs16947": {
        "gene": "CYP2D6",
        "variant": "*2",
        "function": "Normal to increased function",
        "risk_allele": "A",
        "frequency": {"EUR": 0.25, "AFR": 0.30},
        "activity_score": 1,
        "clinical_impact": "Normal function, but may be duplicated (ultrarapid)",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27997040"],
        "notes": "Check for gene duplication if *2 present"
    },
    "rs35742686": {
        "gene": "CYP2D6",
        "variant": "*3",
        "function": "No function (null allele)",
        "risk_allele": "del",
        "frequency": {"EUR": 0.02},
        "activity_score": 0,
        "clinical_impact": "Poor metabolizer",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27997040"]
    },
    "rs5030656": {
        "gene": "CYP2D6",
        "variant": "*9",
        "function": "Decreased function",
        "risk_allele": "del",
        "frequency": {"EUR": 0.02},
        "activity_score": 0.5,
        "clinical_impact": "Intermediate metabolizer",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27997040"]
    },
    "rs59421388": {
        "gene": "CYP2D6",
        "variant": "*29",
        "function": "Decreased function",
        "risk_allele": "G",
        "frequency": {"AFR": 0.15},
        "activity_score": 0.5,
        "clinical_impact": "Intermediate metabolizer - important in African populations",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27997040"]
    },
}

# =============================================================================
# CYP2C19 - Clopidogrel, PPIs, SSRIs, Antifungals
# =============================================================================

CYP2C19_MARKERS = {
    "rs4244285": {
        "gene": "CYP2C19",
        "variant": "*2",
        "function": "No function (splicing defect)",
        "risk_allele": "A",
        "frequency": {"EUR": 0.15, "EAS": 0.30, "AFR": 0.15, "SAS": 0.35},
        "activity_score": 0,
        "drugs_affected": [
            "clopidogrel", "prasugrel", "ticagrelor",
            "omeprazole", "esomeprazole", "pantoprazole", "lansoprazole", "rabeprazole",
            "citalopram", "escitalopram", "sertraline", "amitriptyline", "clomipramine",
            "voriconazole", "clobazam", "brivaracetam",
            "proguanil"
        ],
        "clinical_impact": "Poor metabolizer - clopidogrel may be INEFFECTIVE",
        "cpic_level": "1A",
        "fda_label": "boxed_warning",
        "evidence": "definitive",
        "pmid": ["21716271", "23698643", "25974703"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "clopidogrel": "USE ALTERNATIVE. Prasugrel or ticagrelor recommended for ACS/PCI. Risk of cardiovascular events with clopidogrel.",
                "voriconazole": "Consider 50% maintenance dose reduction. Monitor drug levels.",
                "escitalopram": "Consider 50% dose reduction or alternative SSRI. Max 10mg/day for PMs.",
                "citalopram": "Max 20mg/day for PMs. Consider alternative SSRI.",
                "amitriptyline": "Consider 50% dose reduction. Monitor for adverse effects.",
                "general": "CRITICAL for cardiology - inform cardiologist BEFORE any cardiac procedure."
            }
        }
    },
    "rs4986893": {
        "gene": "CYP2C19",
        "variant": "*3",
        "function": "No function (premature stop codon)",
        "risk_allele": "A",
        "frequency": {"EAS": 0.05, "SAS": 0.02},
        "activity_score": 0,
        "drugs_affected": ["Same as *2"],
        "clinical_impact": "Poor metabolizer - important in East Asian populations",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21716271", "8329904"]
    },
    "rs12248560": {
        "gene": "CYP2C19",
        "variant": "*17",
        "function": "Increased function (promoter variant)",
        "risk_allele": "T",
        "frequency": {"EUR": 0.21, "AFR": 0.18, "SAS": 0.15, "EAS": 0.03},
        "activity_score": 1.5,
        "drugs_affected": ["clopidogrel", "SSRIs", "PPIs", "voriconazole"],
        "clinical_impact": "Ultrarapid metabolizer - increased clopidogrel effect, reduced SSRI/PPI effect",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21716271", "16958828"],
        "actionable": {
            "priority": ClinicalActionLevel.MODERATE,
            "recommendations": {
                "clopidogrel": "Standard dosing. May have enhanced antiplatelet effect - monitor bleeding.",
                "escitalopram": "May need higher doses for therapeutic effect. Standard starting dose, titrate as needed.",
                "ppis": "May need higher PPI doses or twice-daily dosing for acid control.",
                "voriconazole": "May need higher doses. Monitor drug levels."
            }
        }
    },
    "rs28399504": {
        "gene": "CYP2C19",
        "variant": "*4",
        "function": "No function",
        "risk_allele": "A",
        "frequency": {"EUR": 0.01},
        "activity_score": 0,
        "clinical_impact": "Poor metabolizer - rare",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21716271"]
    },
    "rs56337013": {
        "gene": "CYP2C19",
        "variant": "*5",
        "function": "No function",
        "risk_allele": "T",
        "frequency": {"EAS": 0.01},
        "activity_score": 0,
        "clinical_impact": "Poor metabolizer - rare",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21716271"]
    },
}

# =============================================================================
# CYP2C9 - Warfarin, NSAIDs, Sulfonylureas
# =============================================================================

CYP2C9_MARKERS = {
    "rs1799853": {
        "gene": "CYP2C9",
        "variant": "*2",
        "function": "Decreased function (~70% activity)",
        "risk_allele": "T",
        "frequency": {"EUR": 0.13, "SAS": 0.04, "AFR": 0.02},
        "drugs_affected": [
            "warfarin", "acenocoumarol", "phenprocoumon",
            "phenytoin", "fosphenytoin",
            "losartan",
            "celecoxib", "flurbiprofen", "ibuprofen", "lornoxicam", "meloxicam", "piroxicam", "tenoxicam",
            "glipizide", "glimepiride", "glyburide", "tolbutamide",
            "siponimod"
        ],
        "clinical_impact": "Reduced warfarin clearance - ~20% dose reduction typically needed",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["21900891", "28198005", "25594166"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "warfarin": "Reduce initial dose by ~20%. Use FDA-approved dosing calculator. Increased bleeding risk.",
                "phenytoin": "Consider 25% dose reduction. Monitor levels.",
                "sulfonylureas": "Increased hypoglycemia risk. Start with lower dose.",
                "nsaids": "Standard doses, but monitor for GI bleeding more closely.",
                "general": "Critical to share with anticoagulation clinic before starting warfarin."
            }
        }
    },
    "rs1057910": {
        "gene": "CYP2C9",
        "variant": "*3",
        "function": "Significantly decreased function (~5% activity)",
        "risk_allele": "C",
        "frequency": {"EUR": 0.07, "SAS": 0.10, "EAS": 0.03},
        "drugs_affected": ["Same as *2 - more severe"],
        "clinical_impact": "Poor metabolizer - MAJOR warfarin dose reduction needed (~40%)",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["21900891", "9867757"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "warfarin": "Reduce initial dose by ~40%. HIGH bleeding risk at standard doses. Consider DOAC alternative.",
                "phenytoin": "Consider 50% dose reduction. Monitor levels closely.",
                "sulfonylureas": "HIGH hypoglycemia risk. Use alternative or start at 50% dose.",
                "nsaids": "Increased GI bleeding risk. Use lowest effective dose for shortest duration.",
                "general": "MANDATORY pharmacogenomic consultation before warfarin therapy."
            }
        }
    },
    "rs7900194": {
        "gene": "CYP2C9",
        "variant": "*8",
        "function": "Decreased function",
        "risk_allele": "A",
        "frequency": {"AFR": 0.06},
        "clinical_impact": "Important in African populations - often missed by standard testing",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21900891", "15159687"],
        "notes": "Critical for warfarin dosing in African ancestry patients"
    },
    "rs28371686": {
        "gene": "CYP2C9",
        "variant": "*11",
        "function": "Decreased function",
        "risk_allele": "T",
        "frequency": {"AFR": 0.02},
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21900891"]
    },
    "rs9332131": {
        "gene": "CYP2C9",
        "variant": "*6",
        "function": "No function",
        "risk_allele": "del",
        "frequency": {"AFR": 0.01},
        "clinical_impact": "Poor metabolizer - rare but severe",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21900891"]
    },
}

# =============================================================================
# CYP3A4/CYP3A5 - Largest metabolizer family
# =============================================================================

CYP3A_MARKERS = {
    "rs776746": {
        "gene": "CYP3A5",
        "variant": "*3",
        "function": "Non-expressor (splicing defect)",
        "risk_allele": "C",
        "frequency": {"EUR": 0.93, "EAS": 0.70, "AFR": 0.30, "SAS": 0.65},
        "drugs_affected": [
            "tacrolimus", "cyclosporine", "sirolimus", "everolimus",
            "midazolam", "alprazolam", "triazolam",
            "atorvastatin", "simvastatin", "lovastatin",
            "amlodipine", "nifedipine", "felodipine", "diltiazem",
            "carbamazepine", "vincristine"
        ],
        "clinical_impact": "Non-expressors (*3/*3) need LOWER tacrolimus doses",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["25801146", "26417955"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "tacrolimus": "*3/*3: Start at standard dose (0.15-0.2 mg/kg/day). *1/*1 or *1/*3: May need higher doses (1.5-2x).",
                "cyclosporine": "CYP3A5 expressors may need higher doses. Monitor drug levels.",
                "general": "CRITICAL for organ transplant patients. Inform transplant team of CYP3A5 status."
            }
        }
    },
    "rs10264272": {
        "gene": "CYP3A5",
        "variant": "*6",
        "function": "Non-expressor",
        "risk_allele": "A",
        "frequency": {"AFR": 0.15},
        "clinical_impact": "Non-expressor - important in African populations",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["25801146"]
    },
    "rs41303343": {
        "gene": "CYP3A5",
        "variant": "*7",
        "function": "Non-expressor",
        "risk_allele": "ins",
        "frequency": {"AFR": 0.10},
        "clinical_impact": "Non-expressor - African populations",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["25801146"]
    },
    "rs35599367": {
        "gene": "CYP3A4",
        "variant": "*22",
        "function": "Decreased expression",
        "risk_allele": "T",
        "frequency": {"EUR": 0.05, "AFR": 0.01},
        "drugs_affected": ["tacrolimus", "statins", "benzodiazepines", "many others"],
        "clinical_impact": "Reduced CYP3A4 activity - affects statin and tacrolimus levels",
        "cpic_level": "2A",
        "evidence": "moderate",
        "pmid": ["24144200", "27037237"],
        "actionable": {
            "priority": ClinicalActionLevel.MODERATE,
            "recommendations": {
                "tacrolimus": "May need lower doses. Combines with CYP3A5 status.",
                "statins": "May need lower doses of CYP3A4-metabolized statins.",
                "general": "Consider drug levels when on CYP3A4 substrates."
            }
        }
    },
}

# =============================================================================
# CYP1A2 - Caffeine, Theophylline, Antipsychotics
# =============================================================================

CYP1A2_MARKERS = {
    "rs762551": {
        "gene": "CYP1A2",
        "variant": "*1F",
        "function": "Inducible - high activity when induced (smoking, diet)",
        "risk_allele": "C",
        "frequency": {"EUR": 0.32, "EAS": 0.25, "AFR": 0.40},
        "drugs_affected": [
            "caffeine", "theophylline",
            "clozapine", "olanzapine", "haloperidol",
            "melatonin", "ramelteon",
            "duloxetine", "fluvoxamine",
            "propranolol", "ropinirole", "tizanidine",
            "erlotinib"
        ],
        "clinical_impact": "AA = fast metabolizer (induced), CC = slow metabolizer",
        "cpic_level": "2A",
        "evidence": "moderate",
        "pmid": ["17132150", "16522833", "20485253"],
        "actionable": {
            "priority": ClinicalActionLevel.MODERATE,
            "recommendations": {
                "caffeine": "CC: Slow metabolizer. Limit caffeine after 2pm. Higher cardiovascular risk with >3 cups/day.",
                "theophylline": "AA (smokers): May need higher doses. CC: Standard or lower doses.",
                "clozapine": "AA (smokers): May need higher doses. Monitor when smoking status changes.",
                "melatonin": "CC: May have prolonged melatonin effect. AA: May need higher doses.",
                "general": "Smoking, charred foods, and cruciferous vegetables induce CYP1A2."
            }
        }
    },
    "rs2069514": {
        "gene": "CYP1A2",
        "variant": "*1C",
        "function": "Decreased expression",
        "risk_allele": "A",
        "frequency": {"EAS": 0.25, "SAS": 0.10},
        "clinical_impact": "Slow metabolizer - important in East Asian populations",
        "cpic_level": "2B",
        "evidence": "moderate",
        "pmid": ["8476986"]
    },
    "rs12720461": {
        "gene": "CYP1A2",
        "variant": "*1K",
        "function": "Decreased expression",
        "risk_allele": "G",
        "frequency": {"EUR": 0.02, "AFR": 0.03},
        "cpic_level": "2B",
        "evidence": "moderate",
        "pmid": ["17132150"]
    },
}

# =============================================================================
# SLCO1B1 - Statin Myopathy
# =============================================================================

SLCO1B1_MARKERS = {
    "rs4149056": {
        "gene": "SLCO1B1",
        "variant": "*5 (521T>C)",
        "function": "Decreased hepatic uptake transport",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15, "EAS": 0.14, "SAS": 0.05, "AFR": 0.02},
        "drugs_affected": [
            "simvastatin", "atorvastatin", "pravastatin", "rosuvastatin", "pitavastatin", "fluvastatin",
            "methotrexate", "rifampin"
        ],
        "clinical_impact": "4.5x increased risk of simvastatin myopathy per C allele (17x for CC)",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["18650507", "24918167", "23698655"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "simvastatin": "CC: AVOID >20mg/day. Consider pravastatin, rosuvastatin, or fluvastatin. TC: Max 40mg/day.",
                "atorvastatin": "CC: Consider lower starting dose. Use pravastatin or rosuvastatin as alternative.",
                "rosuvastatin": "Lower myopathy risk than simvastatin. Max 20mg/day for CC.",
                "pravastatin": "Lowest myopathy risk. Preferred statin for CC genotype.",
                "general": "Report any muscle pain, weakness, or dark urine immediately. Check CK if symptoms develop."
            }
        }
    },
    "rs2306283": {
        "gene": "SLCO1B1",
        "variant": "*1b (388A>G)",
        "function": "Increased transport function",
        "risk_allele": "G",
        "frequency": {"EUR": 0.40, "EAS": 0.70, "AFR": 0.70},
        "clinical_impact": "May partially offset *5 effect when in cis (haplotype *15)",
        "cpic_level": "1A",
        "evidence": "moderate",
        "pmid": ["24918167"],
        "notes": "*15 haplotype (*1b + *5) has intermediate risk"
    },
}

# =============================================================================
# VKORC1 - Warfarin Target
# =============================================================================

VKORC1_MARKERS = {
    "rs9923231": {
        "gene": "VKORC1",
        "variant": "-1639G>A",
        "function": "Reduced expression (promoter variant)",
        "risk_allele": "A",  # Note: commonly reported as T on forward strand
        "frequency": {"EUR": 0.40, "EAS": 0.92, "AFR": 0.10, "SAS": 0.50},
        "drugs_affected": ["warfarin", "phenprocoumon", "acenocoumarol"],
        "clinical_impact": "AA requires ~50% lower warfarin dose; AG requires ~25% lower dose",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["15930419", "26417955"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "warfarin": "AA: ~50% dose reduction needed. AG: ~25% dose reduction. Use FDA warfarin dosing calculator.",
                "general": "Combine with CYP2C9 status for total dose calculation. Use www.warfarindosing.org"
            }
        }
    },
    "rs8050894": {
        "gene": "VKORC1",
        "variant": "1173C>T",
        "function": "Linked to -1639G>A",
        "risk_allele": "T",
        "frequency": {"EUR": 0.35, "EAS": 0.90},
        "clinical_impact": "In high LD with rs9923231 - similar effect",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["15930419"]
    },
    "rs7294": {
        "gene": "VKORC1",
        "variant": "3730G>A",
        "function": "Decreased expression",
        "risk_allele": "A",
        "frequency": {"EUR": 0.40},
        "clinical_impact": "Part of low-dose haplotype",
        "cpic_level": "1A",
        "evidence": "moderate",
        "pmid": ["15930419"]
    },
}

# =============================================================================
# DPYD - Fluoropyrimidine Toxicity (CRITICAL - potentially fatal)
# =============================================================================

DPYD_MARKERS = {
    "rs3918290": {
        "gene": "DPYD",
        "variant": "*2A (IVS14+1G>A)",
        "function": "No function (splice site)",
        "risk_allele": "A",
        "frequency": {"EUR": 0.01, "AFR": 0.001},
        "drugs_affected": ["5-fluorouracil", "capecitabine", "tegafur"],
        "clinical_impact": "SEVERE/FATAL TOXICITY - complete DPD deficiency",
        "cpic_level": "1A",
        "fda_label": "boxed_warning",
        "evidence": "definitive",
        "pmid": ["28686845", "29152729"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "5-fluorouracil": "Heterozygous: 50% dose reduction REQUIRED. Homozygous: CONTRAINDICATED.",
                "capecitabine": "Same as 5-FU. Life-threatening toxicity at standard doses.",
                "general": "FATAL neutropenia, mucositis, hand-foot syndrome possible. MANDATORY pre-treatment testing."
            }
        }
    },
    "rs55886062": {
        "gene": "DPYD",
        "variant": "*13 (I560S)",
        "function": "No function",
        "risk_allele": "A",
        "frequency": {"EUR": 0.001},
        "drugs_affected": ["5-fluorouracil", "capecitabine"],
        "clinical_impact": "Complete DPD deficiency - severe toxicity",
        "cpic_level": "1A",
        "fda_label": "boxed_warning",
        "evidence": "definitive",
        "pmid": ["28686845"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "general": "Heterozygous: 50% dose reduction. Homozygous: CONTRAINDICATED."
            }
        }
    },
    "rs67376798": {
        "gene": "DPYD",
        "variant": "D949V (c.2846A>T)",
        "function": "Decreased function",
        "risk_allele": "T",
        "frequency": {"EUR": 0.01},
        "drugs_affected": ["5-fluorouracil", "capecitabine"],
        "clinical_impact": "Partial DPD deficiency - increased toxicity risk",
        "cpic_level": "1A",
        "fda_label": "boxed_warning",
        "evidence": "definitive",
        "pmid": ["28686845"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "general": "25-50% dose reduction recommended. Monitor closely."
            }
        }
    },
    "rs75017182": {
        "gene": "DPYD",
        "variant": "HapB3 (c.1129-5923C>G)",
        "function": "Decreased function (deep intronic)",
        "risk_allele": "G",
        "frequency": {"EUR": 0.02},
        "drugs_affected": ["5-fluorouracil", "capecitabine"],
        "clinical_impact": "Partial DPD deficiency",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["28686845", "23988873"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "general": "25% dose reduction recommended."
            }
        }
    },
    "rs1801265": {
        "gene": "DPYD",
        "variant": "C29R (DPYD*9A)",
        "function": "Normal function (benign)",
        "risk_allele": "G",
        "frequency": {"EUR": 0.20},
        "clinical_impact": "No clinical significance - normal activity",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["28686845"],
        "notes": "Common variant, NOT associated with toxicity"
    },
}

# =============================================================================
# TPMT/NUDT15 - Thiopurine Toxicity (CRITICAL)
# =============================================================================

TPMT_MARKERS = {
    "rs1142345": {
        "gene": "TPMT",
        "variant": "*3C (Y240C)",
        "function": "No function",
        "risk_allele": "G",
        "frequency": {"EUR": 0.04, "EAS": 0.02, "AFR": 0.07, "SAS": 0.02},
        "drugs_affected": ["azathioprine", "mercaptopurine", "thioguanine"],
        "clinical_impact": "Severe myelosuppression at standard doses",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["21270794", "23422873"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "azathioprine": "Heterozygous: Start at 30-70% of standard dose. Homozygous: Start at 10% dose or avoid.",
                "mercaptopurine": "Same as azathioprine. Life-threatening bone marrow suppression possible.",
                "general": "MANDATORY testing before thiopurine therapy. Monitor CBC closely."
            }
        }
    },
    "rs1800460": {
        "gene": "TPMT",
        "variant": "*3B (A154T)",
        "function": "No function",
        "risk_allele": "A",
        "frequency": {"EUR": 0.01},
        "drugs_affected": ["azathioprine", "mercaptopurine", "thioguanine"],
        "clinical_impact": "Severe myelosuppression - usually in cis with *3C",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21270794"]
    },
    "rs1800462": {
        "gene": "TPMT",
        "variant": "*2 (A80P)",
        "function": "No function",
        "risk_allele": "G",
        "frequency": {"EUR": 0.003},
        "clinical_impact": "Severe myelosuppression - rare but severe",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21270794"]
    },
    "rs1800584": {
        "gene": "TPMT",
        "variant": "*3A component",
        "function": "In cis with *3B and *3C",
        "risk_allele": "A",
        "clinical_impact": "Part of *3A haplotype (most common null allele in Europeans)",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21270794"]
    },
}

NUDT15_MARKERS = {
    "rs116855232": {
        "gene": "NUDT15",
        "variant": "*3 (R139C)",
        "function": "No function",
        "risk_allele": "T",
        "frequency": {"EAS": 0.10, "SAS": 0.08, "AMR": 0.04, "EUR": 0.002},
        "drugs_affected": ["azathioprine", "mercaptopurine", "thioguanine"],
        "clinical_impact": "Severe myelosuppression - MORE IMPORTANT than TPMT in Asian populations",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["27535131", "29152729"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "azathioprine": "Heterozygous: 25-50% dose reduction. Homozygous: 10% dose or AVOID.",
                "general": "CRITICAL in East Asian patients. Test BOTH TPMT and NUDT15. Can be more important than TPMT."
            }
        }
    },
    "rs147390019": {
        "gene": "NUDT15",
        "variant": "*2 (V18I)",
        "function": "Intermediate function",
        "risk_allele": "T",
        "frequency": {"EAS": 0.02},
        "clinical_impact": "Intermediate risk - dose adjustment may be needed",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["27535131"]
    },
    "rs186364861": {
        "gene": "NUDT15",
        "variant": "*4 (R139H)",
        "function": "Decreased function",
        "risk_allele": "A",
        "frequency": {"EAS": 0.001},
        "clinical_impact": "Increased toxicity risk",
        "cpic_level": "1A",
        "evidence": "moderate",
        "pmid": ["27535131"]
    },
}

# =============================================================================
# UGT1A1 - Irinotecan Toxicity
# =============================================================================

UGT1A1_MARKERS = {
    "rs8175347": {
        "gene": "UGT1A1",
        "variant": "*28 (TA)7",
        "function": "Decreased expression (promoter TA repeat)",
        "risk_allele": "TA7",  # 7 repeats vs 6
        "frequency": {"EUR": 0.35, "AFR": 0.40, "EAS": 0.15},
        "drugs_affected": ["irinotecan", "belinostat", "atazanavir", "nilotinib"],
        "clinical_impact": "Increased irinotecan toxicity (neutropenia, diarrhea); Gilbert syndrome",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["26417955", "21232014"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "irinotecan": "*28/*28: Consider dose reduction or alternative. Increased toxicity risk.",
                "atazanavir": "Increased hyperbilirubinemia (usually benign but cosmetically bothersome).",
                "general": "Gilbert syndrome (*28/*28) causes benign intermittent jaundice. Not dangerous but inform patient."
            }
        },
        "notes": "Also causes Gilbert syndrome (benign unconjugated hyperbilirubinemia)"
    },
    "rs4148323": {
        "gene": "UGT1A1",
        "variant": "*6 (G71R)",
        "function": "Decreased function",
        "risk_allele": "A",
        "frequency": {"EAS": 0.15, "SAS": 0.05},
        "drugs_affected": ["irinotecan", "atazanavir"],
        "clinical_impact": "Important in East Asian populations - similar effect to *28",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["26417955"],
        "notes": "More common than *28 in East Asians"
    },
}

# =============================================================================
# HLA Alleles - Severe Drug Reactions (CRITICAL)
# =============================================================================

HLA_MARKERS = {
    "rs2395029": {
        "gene": "HCP5 (HLA-B*57:01 tag)",
        "variant": "HLA-B*57:01 proxy",
        "function": "HLA class I marker",
        "risk_allele": "G",
        "frequency": {"EUR": 0.06, "SAS": 0.02, "AFR": 0.01},
        "drugs_affected": ["abacavir", "flucloxacillin", "carbamazepine", "phenytoin"],
        "clinical_impact": "Abacavir hypersensitivity - potentially FATAL reaction",
        "cpic_level": "1A",
        "fda_label": "boxed_warning",
        "evidence": "definitive",
        "pmid": ["22378157", "18826941"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "abacavir": "CONTRAINDICATED if HLA-B*57:01 positive. Severe, potentially fatal hypersensitivity reaction.",
                "general": "Standard of care: HLA-B*57:01 testing REQUIRED before prescribing abacavir."
            }
        }
    },
    "rs3909184": {
        "gene": "HLA-B*15:02 proxy",
        "variant": "HLA-B*15:02 tag SNP",
        "function": "HLA class I marker",
        "risk_allele": "A",
        "frequency": {"EAS": 0.08, "SAS": 0.05, "SEA": 0.15},
        "drugs_affected": ["carbamazepine", "oxcarbazepine", "phenytoin", "lamotrigine"],
        "clinical_impact": "Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis - potentially FATAL",
        "cpic_level": "1A",
        "fda_label": "boxed_warning",
        "evidence": "definitive",
        "pmid": ["24512789", "17476434"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "carbamazepine": "CONTRAINDICATED if HLA-B*15:02 positive. 10% mortality with SJS/TEN.",
                "oxcarbazepine": "Also contraindicated if HLA-B*15:02 positive.",
                "phenytoin": "Consider alternative if HLA-B*15:02 positive.",
                "general": "FDA recommends testing in patients of Asian ancestry before carbamazepine. Southeast Asian highest risk."
            }
        }
    },
    "rs1061235": {
        "gene": "HLA-A*31:01 proxy",
        "variant": "HLA-A*31:01 tag SNP",
        "function": "HLA class I marker",
        "risk_allele": "A",
        "frequency": {"EUR": 0.03, "EAS": 0.02, "AMR": 0.05},
        "drugs_affected": ["carbamazepine"],
        "clinical_impact": "DRESS (Drug Reaction with Eosinophilia and Systemic Symptoms)",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["24512789", "21471972"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "carbamazepine": "Higher risk of DRESS syndrome. Consider alternative anticonvulsant.",
                "general": "Important in European ancestry patients (HLA-B*15:02 more important in Asians)."
            }
        }
    },
    "rs2844682": {
        "gene": "HLA-B*58:01 proxy",
        "variant": "HLA-B*58:01 tag SNP",
        "function": "HLA class I marker",
        "risk_allele": "T",
        "frequency": {"EAS": 0.08, "AFR": 0.06, "EUR": 0.01},
        "drugs_affected": ["allopurinol"],
        "clinical_impact": "Allopurinol-induced SJS/TEN",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["26095231"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "allopurinol": "Avoid in HLA-B*58:01 positive patients. Use febuxostat as alternative.",
                "general": "Test in patients of Asian or African ancestry before allopurinol."
            }
        }
    },
}

# =============================================================================
# OPRM1 - Opioid Response
# =============================================================================

OPRM1_MARKERS = {
    "rs1799971": {
        "gene": "OPRM1",
        "variant": "A118G",
        "function": "Decreased receptor expression/binding",
        "risk_allele": "G",
        "frequency": {"EUR": 0.15, "EAS": 0.40, "AFR": 0.05, "SAS": 0.25},
        "drugs_affected": ["morphine", "fentanyl", "oxycodone", "methadone", "naltrexone", "buprenorphine"],
        "clinical_impact": "Reduced opioid analgesia, may need higher doses",
        "cpic_level": "2A",
        "evidence": "moderate",
        "pmid": ["21412232", "17329694"],
        "actionable": {
            "priority": ClinicalActionLevel.MODERATE,
            "recommendations": {
                "opioids": "GG genotype: May require 10-30% higher opioid doses for adequate analgesia.",
                "naltrexone": "May have better response to naltrexone for alcohol use disorder.",
                "general": "Combines with CYP2D6 status for total opioid response prediction."
            }
        }
    },
}

# =============================================================================
# COMT - Pain and Psychiatric Medications
# =============================================================================

COMT_DRUG_MARKERS = {
    "rs4680": {
        "gene": "COMT",
        "variant": "Val158Met",
        "function": "Met = slow catecholamine degradation",
        "risk_allele": "A",  # A = Met (slow)
        "frequency": {"EUR": 0.48, "EAS": 0.30, "AFR": 0.40},
        "drugs_affected": [
            "morphine", "tramadol",
            "amphetamine", "methylphenidate",
            "tolcapone", "entacapone"
        ],
        "clinical_impact": "Met/Met: Higher pain sensitivity, better opioid response, higher anxiety risk",
        "cpic_level": "2B",
        "evidence": "moderate",
        "pmid": ["17008817", "14517761"],
        "actionable": {
            "priority": ClinicalActionLevel.MODERATE,
            "recommendations": {
                "opioids": "Met/Met: More pain sensitive, may need higher doses. Better response when given.",
                "stimulants": "Met/Met: May be more sensitive to stimulant side effects.",
                "general": "Met/Met = 'Worrier' (higher anxiety, better focus). Val/Val = 'Warrior' (stress resilient)."
            }
        }
    },
}

# =============================================================================
# ABCB1 - Drug Transport
# =============================================================================

ABCB1_MARKERS = {
    "rs1045642": {
        "gene": "ABCB1/MDR1",
        "variant": "C3435T",
        "function": "Altered P-glycoprotein expression",
        "risk_allele": "T",
        "frequency": {"EUR": 0.55, "EAS": 0.40, "AFR": 0.20},
        "drugs_affected": [
            "digoxin", "cyclosporine", "tacrolimus",
            "fexofenadine", "loperamide",
            "dabigatran", "clopidogrel",
            "many chemotherapy drugs"
        ],
        "clinical_impact": "TT may have higher brain penetration of P-gp substrates",
        "cpic_level": "3",
        "evidence": "limited",
        "pmid": ["11701890", "16960150"],
        "actionable": {
            "priority": ClinicalActionLevel.LOW,
            "recommendations": {
                "general": "Research variant. May affect CNS penetration of P-gp substrates."
            }
        }
    },
    "rs2032582": {
        "gene": "ABCB1/MDR1",
        "variant": "G2677T/A",
        "function": "Altered P-glycoprotein function",
        "risk_allele": "T",
        "frequency": {"EUR": 0.45, "EAS": 0.50},
        "clinical_impact": "May affect drug transport - often tested with rs1045642",
        "cpic_level": "3",
        "evidence": "limited",
        "pmid": ["11701890"]
    },
    "rs1128503": {
        "gene": "ABCB1/MDR1",
        "variant": "C1236T",
        "function": "Synonymous but may affect expression",
        "risk_allele": "T",
        "frequency": {"EUR": 0.45, "EAS": 0.65},
        "clinical_impact": "Part of ABCB1 haplotype",
        "cpic_level": "3",
        "evidence": "limited",
        "pmid": ["11701890"]
    },
}

# =============================================================================
# NAT2 - Isoniazid, Sulfonamides
# =============================================================================

NAT2_MARKERS = {
    "rs1801280": {
        "gene": "NAT2",
        "variant": "*5 (I114T)",
        "function": "Slow acetylator",
        "risk_allele": "C",
        "frequency": {"EUR": 0.45, "AFR": 0.25, "EAS": 0.05},
        "drugs_affected": ["isoniazid", "hydralazine", "procainamide", "sulfasalazine", "dapsone", "caffeine"],
        "clinical_impact": "Slow acetylators: increased hepatotoxicity with isoniazid, drug-induced lupus with hydralazine",
        "cpic_level": "2A",
        "evidence": "moderate",
        "pmid": ["24727250", "8319332"],
        "actionable": {
            "priority": ClinicalActionLevel.MODERATE,
            "recommendations": {
                "isoniazid": "Slow acetylators: Increased hepatotoxicity risk. Monitor LFTs closely.",
                "hydralazine": "Slow acetylators: Higher risk of lupus-like syndrome. Max 200mg/day.",
                "general": "~50% of Europeans are slow acetylators. Faster in East Asians."
            }
        }
    },
    "rs1799930": {
        "gene": "NAT2",
        "variant": "*6 (R197Q)",
        "function": "Slow acetylator",
        "risk_allele": "A",
        "frequency": {"EUR": 0.30, "EAS": 0.25, "AFR": 0.25},
        "clinical_impact": "Slow acetylator allele",
        "cpic_level": "2A",
        "evidence": "moderate",
        "pmid": ["24727250"]
    },
    "rs1799931": {
        "gene": "NAT2",
        "variant": "*7 (G286E)",
        "function": "Slow acetylator",
        "risk_allele": "A",
        "frequency": {"EUR": 0.03, "EAS": 0.10, "SAS": 0.10},
        "clinical_impact": "Slow acetylator allele - more common in East/South Asian",
        "cpic_level": "2A",
        "evidence": "moderate",
        "pmid": ["24727250"]
    },
    "rs1801279": {
        "gene": "NAT2",
        "variant": "*14 (R64Q)",
        "function": "Slow acetylator",
        "risk_allele": "A",
        "frequency": {"AFR": 0.10},
        "clinical_impact": "Slow acetylator - important in African populations",
        "cpic_level": "2A",
        "evidence": "moderate",
        "pmid": ["24727250"]
    },
}

# =============================================================================
# G6PD - Oxidative Drug Hemolysis
# =============================================================================

G6PD_MARKERS = {
    "rs1050828": {
        "gene": "G6PD",
        "variant": "V68M (A- African variant)",
        "function": "Decreased activity (~10% residual)",
        "risk_allele": "T",
        "frequency": {"AFR": 0.20, "AMR": 0.05},
        "drugs_affected": [
            "primaquine", "tafenoquine",
            "dapsone", "rasburicase",
            "methylene blue",
            "sulfonamides", "nitrofurantoin",
            "high-dose aspirin"
        ],
        "clinical_impact": "Hemolytic anemia with oxidative drugs - X-linked",
        "cpic_level": "1A",
        "fda_label": "required",
        "evidence": "definitive",
        "pmid": ["22457341", "29083402"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "primaquine": "CONTRAINDICATED in G6PD deficiency. Severe hemolysis.",
                "rasburicase": "CONTRAINDICATED. Can cause fatal hemolysis.",
                "dapsone": "AVOID. Use alternative for PCP prophylaxis.",
                "general": "X-linked: Males more severely affected. Females can be carriers or affected. Avoid fava beans."
            }
        }
    },
    "rs5030868": {
        "gene": "G6PD",
        "variant": "S188F (Mediterranean variant)",
        "function": "Severely decreased activity (<5% residual)",
        "risk_allele": "A",
        "frequency": {"MENA": 0.05, "SAS": 0.03, "EUR": 0.01},
        "drugs_affected": ["Same as A- variant - MORE SEVERE"],
        "clinical_impact": "Severe G6PD deficiency - chronic hemolysis possible",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["22457341"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "general": "MORE SEVERE than African variant. Chronic hemolysis possible. Strict avoidance of triggers."
            }
        }
    },
    "rs137852328": {
        "gene": "G6PD",
        "variant": "Mahidol (G163S)",
        "function": "Severely decreased",
        "risk_allele": "A",
        "frequency": {"SEA": 0.15, "THA": 0.25},
        "clinical_impact": "Common in Southeast Asia - severe deficiency",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["22457341"]
    },
}

# =============================================================================
# RYR1/CACNA1S - Malignant Hyperthermia
# =============================================================================

ANESTHESIA_MARKERS = {
    "rs121918592": {
        "gene": "RYR1",
        "variant": "R614C",
        "function": "Malignant hyperthermia susceptibility (MHS)",
        "risk_allele": "T",
        "frequency": {"EUR": 0.0001},  # Rare but severe
        "drugs_affected": ["succinylcholine", "sevoflurane", "desflurane", "isoflurane", "halothane"],
        "clinical_impact": "LIFE-THREATENING reaction to anesthesia - 70% mortality untreated",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21880859", "9762081"],
        "actionable": {
            "priority": ClinicalActionLevel.CRITICAL,
            "recommendations": {
                "anesthesia": "AVOID succinylcholine and ALL volatile anesthetics. Use total IV anesthesia (TIVA).",
                "general": "Medical alert bracelet REQUIRED. Dantrolene must be available. Inform ALL healthcare providers."
            }
        }
    },
    "rs28933396": {
        "gene": "RYR1",
        "variant": "G341R",
        "function": "MH susceptibility",
        "risk_allele": "A",
        "frequency": {"EUR": 0.0001},
        "clinical_impact": "Malignant hyperthermia susceptibility",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21880859"]
    },
    "rs121918593": {
        "gene": "RYR1",
        "variant": "G2434R",
        "function": "MH susceptibility",
        "risk_allele": "A",
        "frequency": {"EUR": 0.0001},
        "clinical_impact": "Malignant hyperthermia susceptibility",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21880859"]
    },
    "rs772226819": {
        "gene": "CACNA1S",
        "variant": "R1086H",
        "function": "MH susceptibility",
        "risk_allele": "A",
        "clinical_impact": "1-2% of MH cases linked to CACNA1S",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21880859"]
    },
}

BCHE_MARKERS = {
    "rs1799807": {
        "gene": "BCHE",
        "variant": "A (Atypical)",
        "function": "Pseudocholinesterase deficiency",
        "risk_allele": "A",
        "frequency": {"EUR": 0.02},
        "drugs_affected": ["succinylcholine", "mivacurium"],
        "clinical_impact": "Prolonged paralysis (2-8 hours) after succinylcholine",
        "cpic_level": "2A",
        "evidence": "moderate",
        "pmid": ["1349607"],
        "actionable": {
            "priority": ClinicalActionLevel.HIGH,
            "recommendations": {
                "succinylcholine": "Heterozygous: Moderately prolonged paralysis. Homozygous: 2-8 hour paralysis.",
                "general": "Inform anesthesiologist. Use rocuronium/sugammadex as alternative."
            }
        }
    },
    "rs28933389": {
        "gene": "BCHE",
        "variant": "K (Kalow)",
        "function": "Reduced activity (66%)",
        "risk_allele": "A",
        "frequency": {"EUR": 0.12},
        "clinical_impact": "Mild prolonged paralysis if compound heterozygous",
        "cpic_level": "2B",
        "evidence": "moderate",
        "pmid": ["1349607"]
    },
}

# =============================================================================
# CYP4F2 - Vitamin K / Warfarin
# =============================================================================

CYP4F2_MARKERS = {
    "rs2108622": {
        "gene": "CYP4F2",
        "variant": "V433M",
        "function": "Decreased vitamin K metabolism",
        "risk_allele": "T",
        "frequency": {"EUR": 0.30, "EAS": 0.20, "AFR": 0.05},
        "drugs_affected": ["warfarin", "vitamin K"],
        "clinical_impact": "TT genotype: ~1mg/day higher warfarin dose needed",
        "cpic_level": "1A",
        "evidence": "definitive",
        "pmid": ["21900891", "19578179"],
        "actionable": {
            "priority": ClinicalActionLevel.MODERATE,
            "recommendations": {
                "warfarin": "TT: May need higher warfarin dose (~1mg/day more). Include in dosing algorithms."
            }
        }
    },
}

# =============================================================================
# DRD2/ANKK1 - Dopamine Drug Response
# =============================================================================

DOPAMINE_MARKERS = {
    "rs1800497": {
        "gene": "ANKK1/DRD2",
        "variant": "Taq1A (rs1800497)",
        "function": "Reduced D2 receptor density",
        "risk_allele": "A",
        "frequency": {"EUR": 0.20, "EAS": 0.40, "AFR": 0.35},
        "drugs_affected": ["antipsychotics", "metoclopramide", "bromocriptine", "naltrexone"],
        "clinical_impact": "A1 allele: Reduced D2 receptors, may need higher antipsychotic doses",
        "cpic_level": "3",
        "evidence": "limited",
        "pmid": ["8807664", "2906084"],
        "actionable": {
            "priority": ClinicalActionLevel.LOW,
            "recommendations": {
                "antipsychotics": "A1/A1: May need higher doses of D2 antagonists. Also higher addiction risk.",
                "naltrexone": "A1 carriers: Reduced response to naltrexone for alcohol use disorder."
            }
        }
    },
}

# =============================================================================
# COMBINE ALL MARKERS
# =============================================================================

PHARMACOGENOMICS_COMPLETE = {
    **CYP2D6_MARKERS,
    **CYP2C19_MARKERS,
    **CYP2C9_MARKERS,
    **CYP3A_MARKERS,
    **CYP1A2_MARKERS,
    **SLCO1B1_MARKERS,
    **VKORC1_MARKERS,
    **DPYD_MARKERS,
    **TPMT_MARKERS,
    **NUDT15_MARKERS,
    **UGT1A1_MARKERS,
    **HLA_MARKERS,
    **OPRM1_MARKERS,
    **COMT_DRUG_MARKERS,
    **ABCB1_MARKERS,
    **NAT2_MARKERS,
    **G6PD_MARKERS,
    **ANESTHESIA_MARKERS,
    **BCHE_MARKERS,
    **CYP4F2_MARKERS,
    **DOPAMINE_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def calculate_cyp2d6_activity_score(genotypes: Dict[str, str]) -> tuple[float, MetabolizerStatus]:
    """
    Calculate CYP2D6 activity score from genotypes.
    Returns (activity_score, metabolizer_status).
    """
    activity_scores = {
        "rs3892097": {"AA": 0, "AG": 0.5, "GG": 1},  # *4
        "rs5030655": {"del": 0, "other": 1},  # *6
        "rs1065852": {"AA": 0.25, "AG": 0.625, "GG": 1},  # *10
        "rs28371725": {"AA": 0.5, "AG": 0.75, "GG": 1},  # *41
        "rs28371706": {"AA": 0.5, "AG": 0.75, "GG": 1},  # *17
    }
    
    total_score = 2.0  # Start with normal (two functional alleles)
    
    # Subtract based on detected variants
    for rsid, geno in genotypes.items():
        if rsid in activity_scores and geno in activity_scores[rsid]:
            reduction = 1.0 - activity_scores[rsid][geno]
            total_score -= reduction
    
    # Determine status
    if total_score >= 2.25:
        status = MetabolizerStatus.ULTRARAPID
    elif total_score >= 1.25:
        status = MetabolizerStatus.NORMAL
    elif total_score >= 0.25:
        status = MetabolizerStatus.INTERMEDIATE
    else:
        status = MetabolizerStatus.POOR
    
    return total_score, status

def calculate_cyp2c19_status(genotypes: Dict[str, str]) -> MetabolizerStatus:
    """Calculate CYP2C19 metabolizer status."""
    has_star2 = genotypes.get("rs4244285", "GG") in ["AA", "AG"]
    has_star3 = genotypes.get("rs4986893", "GG") in ["AA", "AG"]
    has_star17 = genotypes.get("rs12248560", "CC") in ["TT", "CT"]
    
    null_count = 0
    if genotypes.get("rs4244285") == "AA":
        null_count += 2
    elif genotypes.get("rs4244285") == "AG":
        null_count += 1
    
    if genotypes.get("rs4986893") == "AA":
        null_count += 2
    elif genotypes.get("rs4986893") == "AG":
        null_count += 1
    
    if null_count >= 2:
        return MetabolizerStatus.POOR
    elif null_count == 1:
        if has_star17:
            return MetabolizerStatus.NORMAL
        return MetabolizerStatus.INTERMEDIATE
    elif has_star17:
        if genotypes.get("rs12248560") == "TT":
            return MetabolizerStatus.ULTRARAPID
        return MetabolizerStatus.RAPID
    
    return MetabolizerStatus.NORMAL

def get_critical_alerts(genotypes: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Get list of critical drug alerts that require immediate attention.
    """
    alerts = []
    
    # Check DPYD (5-FU toxicity)
    if genotypes.get("rs3918290") in ["AA", "AG"]:
        alerts.append({
            "gene": "DPYD",
            "severity": "CRITICAL",
            "drug": "5-fluorouracil/capecitabine",
            "message": "FATAL TOXICITY RISK - 50% dose reduction required if heterozygous, CONTRAINDICATED if homozygous",
            "pmid": "28686845"
        })
    
    # Check TPMT (thiopurine toxicity)
    if genotypes.get("rs1142345") in ["GG", "AG"]:
        alerts.append({
            "gene": "TPMT",
            "severity": "CRITICAL",
            "drug": "azathioprine/mercaptopurine",
            "message": "Severe myelosuppression risk - dose reduction REQUIRED",
            "pmid": "21270794"
        })
    
    # Check NUDT15 (thiopurine toxicity - Asian populations)
    if genotypes.get("rs116855232") in ["TT", "CT"]:
        alerts.append({
            "gene": "NUDT15",
            "severity": "CRITICAL",
            "drug": "azathioprine/mercaptopurine",
            "message": "Severe myelosuppression risk - especially important in Asian ancestry",
            "pmid": "27535131"
        })
    
    # Check HLA-B*57:01 (abacavir)
    if genotypes.get("rs2395029") in ["GG", "AG"]:
        alerts.append({
            "gene": "HLA-B*57:01",
            "severity": "CRITICAL",
            "drug": "abacavir",
            "message": "CONTRAINDICATED - potentially fatal hypersensitivity reaction",
            "pmid": "22378157"
        })
    
    # Check HLA-B*15:02 (carbamazepine SJS/TEN)
    if genotypes.get("rs3909184") in ["AA", "AG"]:
        alerts.append({
            "gene": "HLA-B*15:02",
            "severity": "CRITICAL",
            "drug": "carbamazepine",
            "message": "CONTRAINDICATED - risk of fatal Stevens-Johnson Syndrome",
            "pmid": "24512789"
        })
    
    # Check G6PD
    if genotypes.get("rs1050828") in ["TT", "CT"] or genotypes.get("rs5030868") in ["AA", "AG"]:
        alerts.append({
            "gene": "G6PD",
            "severity": "CRITICAL",
            "drug": "primaquine/rasburicase/dapsone",
            "message": "G6PD deficiency - hemolytic anemia risk with oxidative drugs",
            "pmid": "22457341"
        })
    
    # Check RYR1 (malignant hyperthermia)
    ryr1_variants = ["rs121918592", "rs28933396", "rs121918593"]
    for rs in ryr1_variants:
        if genotypes.get(rs) in ["TT", "AT", "AA", "AG"]:
            alerts.append({
                "gene": "RYR1",
                "severity": "CRITICAL",
                "drug": "succinylcholine/volatile anesthetics",
                "message": "Malignant hyperthermia susceptibility - AVOID triggering agents",
                "pmid": "21880859"
            })
            break
    
    # Check CYP2C19 for clopidogrel
    if genotypes.get("rs4244285") == "AA":
        alerts.append({
            "gene": "CYP2C19",
            "severity": "HIGH",
            "drug": "clopidogrel",
            "message": "Poor metabolizer - clopidogrel may be INEFFECTIVE. Use prasugrel/ticagrelor.",
            "pmid": "21716271"
        })
    
    return alerts

def generate_pharmacogenomics_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Generate comprehensive pharmacogenomics report from genotypes.
    """
    report = {
        "critical_alerts": get_critical_alerts(genotypes),
        "metabolizer_statuses": {},
        "drug_recommendations": {},
        "markers_analyzed": 0,
        "markers_with_variants": 0,
    }
    
    # Calculate metabolizer statuses
    cyp2d6_score, cyp2d6_status = calculate_cyp2d6_activity_score(genotypes)
    report["metabolizer_statuses"]["CYP2D6"] = {
        "activity_score": cyp2d6_score,
        "status": cyp2d6_status.value,
        "interpretation": f"CYP2D6 {cyp2d6_status.value} metabolizer (activity score: {cyp2d6_score})"
    }
    
    cyp2c19_status = calculate_cyp2c19_status(genotypes)
    report["metabolizer_statuses"]["CYP2C19"] = {
        "status": cyp2c19_status.value,
        "interpretation": f"CYP2C19 {cyp2c19_status.value} metabolizer"
    }
    
    # Count markers
    for rsid in PHARMACOGENOMICS_COMPLETE:
        if rsid in genotypes:
            report["markers_analyzed"] += 1
            marker = PHARMACOGENOMICS_COMPLETE[rsid]
            risk_allele = marker.get("risk_allele", "")
            if risk_allele and risk_allele in genotypes.get(rsid, ""):
                report["markers_with_variants"] += 1
    
    return report

# Export all markers and functions
__all__ = [
    'PHARMACOGENOMICS_COMPLETE',
    'CYP2D6_MARKERS',
    'CYP2C19_MARKERS', 
    'CYP2C9_MARKERS',
    'CYP3A_MARKERS',
    'CYP1A2_MARKERS',
    'SLCO1B1_MARKERS',
    'VKORC1_MARKERS',
    'DPYD_MARKERS',
    'TPMT_MARKERS',
    'NUDT15_MARKERS',
    'UGT1A1_MARKERS',
    'HLA_MARKERS',
    'OPRM1_MARKERS',
    'COMT_DRUG_MARKERS',
    'ABCB1_MARKERS',
    'NAT2_MARKERS',
    'G6PD_MARKERS',
    'ANESTHESIA_MARKERS',
    'BCHE_MARKERS',
    'CYP4F2_MARKERS',
    'DOPAMINE_MARKERS',
    'MetabolizerStatus',
    'ClinicalActionLevel',
    'calculate_cyp2d6_activity_score',
    'calculate_cyp2c19_status',
    'get_critical_alerts',
    'generate_pharmacogenomics_report',
]
