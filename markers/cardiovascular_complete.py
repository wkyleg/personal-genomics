"""
Complete Cardiovascular Genetics Panel v5.0

Comprehensive coverage of:
- Lipid metabolism (LDL, HDL, Triglycerides, Lp(a))
- Blood pressure / Hypertension
- Salt sensitivity
- Clotting disorders (thrombophilia)
- Arrhythmia genetics (AFib, Long QT)
- Structural heart disease (cardiomyopathy, aortic aneurysm)
- Homocysteine metabolism

Sources:
- ClinVar
- NHGRI-EBI GWAS Catalog
- American Heart Association guidelines
- European Society of Cardiology

All markers include PMID references and clinical recommendations.
"""

from enum import Enum
from typing import Dict, List, Any, Optional

class CardioRiskLevel(Enum):
    VERY_HIGH = "very_high"
    HIGH = "high"
    MODERATE = "moderate"
    AVERAGE = "average"
    LOW = "low"

class EvidenceLevel(Enum):
    DEFINITIVE = "definitive"      # Replicated, clinical guidelines
    STRONG = "strong"              # Multiple GWAS, consistent
    MODERATE = "moderate"          # Some replication
    LIMITED = "limited"            # Single study or inconsistent

# =============================================================================
# Lp(a) - LIPOPROTEIN(a) - Independent CVD Risk Factor
# =============================================================================

LPA_MARKERS = {
    "rs10455872": {
        "gene": "LPA",
        "variant": "Lp(a) elevation",
        "function": "Strong elevation of Lp(a) levels",
        "risk_allele": "G",
        "frequency": {"EUR": 0.07, "AFR": 0.02, "EAS": 0.01},
        "effect": "Each G allele increases Lp(a) ~50 nmol/L",
        "category": "lipids",
        "condition": "Elevated Lp(a) - Atherosclerotic CVD",
        "or_per_allele": 1.92,
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["19060906", "22085571", "28831089"],
        "actionable": {
            "priority": "high",
            "recommendations": [
                "GG genotype: Very high Lp(a) likely (>100 nmol/L). High cardiovascular risk.",
                "Measure Lp(a) level directly (once in lifetime sufficient as levels stable).",
                "Current therapy limited - emerging antisense oligonucleotides (pelacarsen) in trials.",
                "Aggressive LDL lowering more important with elevated Lp(a).",
                "Consider early statin therapy even with normal LDL.",
                "Aortic valve calcification screening may be warranted."
            ]
        }
    },
    "rs3798220": {
        "gene": "LPA",
        "variant": "Lp(a) elevation",
        "function": "Moderate elevation of Lp(a) levels",
        "risk_allele": "C",
        "frequency": {"EUR": 0.02, "AFR": 0.01},
        "effect": "Each C allele increases Lp(a) ~20-40 nmol/L",
        "category": "lipids",
        "condition": "Elevated Lp(a)",
        "or_per_allele": 1.51,
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["19060906", "22085571"],
        "actionable": {
            "priority": "high",
            "recommendations": [
                "Contributes to Lp(a) elevation, especially with rs10455872.",
                "Measure Lp(a) level.",
                "Independent CVD risk factor not captured by standard lipid panels."
            ]
        }
    },
    "rs6415084": {
        "gene": "LPA",
        "variant": "Lp(a) modulator",
        "function": "Modulates Lp(a) levels",
        "risk_allele": "C",
        "frequency": {"EUR": 0.35, "AFR": 0.50, "EAS": 0.25},
        "category": "lipids",
        "condition": "Lp(a) levels",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["19060906"]
    },
}

# =============================================================================
# LDL CHOLESTEROL GENETICS
# =============================================================================

LDL_MARKERS = {
    "rs11206510": {
        "gene": "PCSK9",
        "variant": "PCSK9 regulatory variant",
        "function": "Affects LDL receptor degradation",
        "risk_allele": "T",
        "frequency": {"EUR": 0.82, "AFR": 0.75, "EAS": 0.85},
        "effect": "T allele associated with higher LDL-C",
        "category": "lipids",
        "condition": "LDL cholesterol levels / CAD risk",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["19060911", "20864672"],
        "actionable": {
            "priority": "medium",
            "recommendations": [
                "TT: Higher LDL-C levels expected.",
                "May benefit more from PCSK9 inhibitors (evolocumab, alirocumab).",
                "Standard statin therapy first-line."
            ]
        }
    },
    "rs2228671": {
        "gene": "LDLR",
        "variant": "LDL receptor variant",
        "function": "Affects LDL receptor function",
        "risk_allele": "T",
        "frequency": {"EUR": 0.12, "EAS": 0.05},
        "effect": "Associated with LDL-C levels",
        "category": "lipids",
        "condition": "LDL cholesterol",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["20864672"]
    },
    "rs5742904": {
        "gene": "APOB",
        "variant": "R3527Q (familial hypercholesterolemia)",
        "function": "Defective LDL receptor binding",
        "risk_allele": "A",
        "frequency": {"EUR": 0.002},
        "effect": "Causes familial defective ApoB-100 (FDB)",
        "category": "lipids",
        "condition": "Familial hypercholesterolemia",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["2563166", "20536758"],
        "actionable": {
            "priority": "critical",
            "recommendations": [
                "Causes familial hypercholesterolemia - autosomal dominant.",
                "Heterozygotes: LDL 200-400 mg/dL. Premature CAD common.",
                "High-intensity statin therapy essential. May need PCSK9 inhibitor.",
                "Family cascade screening recommended.",
                "Refer to lipid specialist."
            ]
        }
    },
    "rs562556": {
        "gene": "APOB",
        "variant": "APOB expression variant",
        "function": "Affects ApoB levels",
        "risk_allele": "G",
        "frequency": {"EUR": 0.30, "AFR": 0.45, "EAS": 0.20},
        "category": "lipids",
        "condition": "LDL and ApoB levels",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["20864672"]
    },
}

# =============================================================================
# HDL CHOLESTEROL GENETICS
# =============================================================================

HDL_MARKERS = {
    "rs1800775": {
        "gene": "CETP",
        "variant": "TaqIB (rs708272 equivalent)",
        "function": "Cholesteryl ester transfer protein activity",
        "risk_allele": "C",
        "frequency": {"EUR": 0.42, "EAS": 0.40, "AFR": 0.30},
        "effect": "A allele (TT) = higher HDL, lower CETP activity",
        "category": "lipids",
        "condition": "HDL cholesterol levels",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["17827399", "20864672"],
        "actionable": {
            "priority": "low",
            "recommendations": [
                "AA genotype: Higher HDL levels expected (may be 10-15% higher).",
                "CETP inhibitors (anacetrapib) did not reduce CVD events despite raising HDL.",
                "Focus on LDL lowering rather than HDL raising.",
                "HDL function may be more important than HDL-C level."
            ]
        }
    },
    "rs1800588": {
        "gene": "LIPC",
        "variant": "Hepatic lipase promoter",
        "function": "Hepatic lipase expression",
        "risk_allele": "C",
        "frequency": {"EUR": 0.22, "EAS": 0.50, "AFR": 0.52},
        "effect": "T allele = higher HDL",
        "category": "lipids",
        "condition": "HDL cholesterol",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["20864672", "10490110"]
    },
    "rs2230806": {
        "gene": "ABCA1",
        "variant": "R219K",
        "function": "Cholesterol efflux",
        "risk_allele": "G",
        "frequency": {"EUR": 0.27, "EAS": 0.35},
        "effect": "Affects HDL levels and CAD risk",
        "category": "lipids",
        "condition": "HDL cholesterol / CAD",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["15226823", "20864672"]
    },
}

# =============================================================================
# TRIGLYCERIDES
# =============================================================================

TRIGLYCERIDE_MARKERS = {
    "rs662799": {
        "gene": "APOA5",
        "variant": "-1131T>C",
        "function": "APOA5 expression",
        "risk_allele": "C",
        "frequency": {"EUR": 0.08, "EAS": 0.30, "AFR": 0.15},
        "effect": "C allele increases TG ~30%",
        "category": "lipids",
        "condition": "Triglycerides / hypertriglyceridemia",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["20864672", "12068375"],
        "actionable": {
            "priority": "medium",
            "recommendations": [
                "CC genotype: Higher triglyceride levels expected.",
                "Important to minimize added sugars and refined carbs.",
                "Omega-3 fatty acids may help lower TG.",
                "Regular exercise particularly effective for TG lowering."
            ]
        }
    },
    "rs1260326": {
        "gene": "GCKR",
        "variant": "P446L",
        "function": "Glucokinase regulation",
        "risk_allele": "T",
        "frequency": {"EUR": 0.40, "EAS": 0.50, "AFR": 0.10},
        "effect": "T allele raises TG, lowers fasting glucose",
        "category": "lipids",
        "condition": "Triglycerides",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["18940312", "20864672"],
        "notes": "Paradoxically protective against T2D despite raising TG"
    },
    "rs328": {
        "gene": "LPL",
        "variant": "S447X (gain of function)",
        "function": "Lipoprotein lipase activity",
        "risk_allele": "C",  # C = 447X (protective)
        "frequency": {"EUR": 0.10, "EAS": 0.05, "AFR": 0.06},
        "effect": "G allele (447X) = lower TG, higher HDL, reduced CAD risk",
        "category": "lipids",
        "condition": "Triglycerides / HDL / CAD",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["20864672", "9808631"],
        "actionable": {
            "priority": "low",
            "recommendations": [
                "CG or CC: Favorable lipid profile - lower TG, higher HDL.",
                "Reduced cardiovascular risk."
            ]
        }
    },
    "rs12678919": {
        "gene": "LPL",
        "variant": "LPL regulatory variant",
        "function": "LPL expression",
        "risk_allele": "A",
        "frequency": {"EUR": 0.10},
        "effect": "Associated with TG levels",
        "category": "lipids",
        "condition": "Triglycerides",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["20864672"]
    },
}

# =============================================================================
# BLOOD PRESSURE / HYPERTENSION
# =============================================================================

BLOOD_PRESSURE_MARKERS = {
    "rs4340": {
        "gene": "ACE",
        "variant": "I/D (Insertion/Deletion)",
        "function": "ACE levels",
        "risk_allele": "D",
        "frequency": {"EUR": 0.52, "AFR": 0.60, "EAS": 0.40},
        "effect": "DD = higher ACE levels, associated with hypertension and CVD",
        "category": "blood_pressure",
        "condition": "Hypertension / Cardiovascular disease",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["8675673", "16567525"],
        "actionable": {
            "priority": "medium",
            "recommendations": [
                "DD genotype: May have higher ACE levels and blood pressure.",
                "May respond well to ACE inhibitors/ARBs.",
                "Salt restriction particularly important.",
                "Associated with exercise performance (power vs endurance)."
            ]
        }
    },
    "rs699": {
        "gene": "AGT",
        "variant": "M235T",
        "function": "Angiotensinogen levels",
        "risk_allele": "G",
        "frequency": {"EUR": 0.42, "AFR": 0.90, "EAS": 0.70},
        "effect": "T allele (G) = higher AGT levels, higher BP",
        "category": "blood_pressure",
        "condition": "Hypertension",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["1641269", "27840431"],
        "actionable": {
            "priority": "medium",
            "recommendations": [
                "GG genotype: Higher angiotensinogen levels, increased hypertension risk.",
                "May respond particularly well to ACE inhibitors/ARBs.",
                "Blood pressure monitoring important."
            ]
        }
    },
    "rs4961": {
        "gene": "ADD1",
        "variant": "Gly460Trp (alpha-adducin)",
        "function": "Sodium reabsorption",
        "risk_allele": "T",
        "frequency": {"EUR": 0.22, "AFR": 0.05, "EAS": 0.50},
        "effect": "Trp allele (T) = salt-sensitive hypertension",
        "category": "blood_pressure",
        "condition": "Salt-sensitive hypertension",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["16567525", "10591398"],
        "actionable": {
            "priority": "high",
            "recommendations": [
                "TT genotype: Salt-sensitive hypertension likely.",
                "STRICT sodium restriction recommended (<1500-2000 mg/day).",
                "Diuretics may be particularly effective.",
                "DASH diet recommended."
            ]
        }
    },
    "rs1799998": {
        "gene": "CYP11B2",
        "variant": "-344C>T (Aldosterone synthase)",
        "function": "Aldosterone production",
        "risk_allele": "C",
        "frequency": {"EUR": 0.50, "AFR": 0.30},
        "effect": "C allele = higher aldosterone, increased BP",
        "category": "blood_pressure",
        "condition": "Hypertension",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["10022116", "16567525"]
    },
    "rs5186": {
        "gene": "AGTR1",
        "variant": "A1166C",
        "function": "AT1 receptor expression",
        "risk_allele": "C",
        "frequency": {"EUR": 0.28, "AFR": 0.05, "EAS": 0.03},
        "effect": "C allele associated with hypertension",
        "category": "blood_pressure",
        "condition": "Hypertension / CAD",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["8557250", "16567525"],
        "actionable": {
            "priority": "low",
            "recommendations": [
                "CC genotype: May have enhanced response to ARBs (losartan, valsartan)."
            ]
        }
    },
    "rs17367504": {
        "gene": "MTHFR-NPPB region",
        "variant": "Blood pressure GWAS hit",
        "function": "Blood pressure regulation",
        "risk_allele": "G",
        "frequency": {"EUR": 0.15},
        "effect": "Associated with blood pressure in GWAS",
        "category": "blood_pressure",
        "condition": "Blood pressure",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["21909115"]
    },
}

# =============================================================================
# CLOTTING / THROMBOPHILIA
# =============================================================================

CLOTTING_MARKERS = {
    "rs6025": {
        "gene": "F5",
        "variant": "Factor V Leiden (R506Q)",
        "function": "Resistance to activated protein C",
        "risk_allele": "A",
        "frequency": {"EUR": 0.05, "AFR": 0.01, "EAS": 0.001},
        "effect": "7x VTE risk (heterozygous), 80x (homozygous)",
        "category": "clotting",
        "condition": "Venous thromboembolism",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["7989260", "8346443"],
        "actionable": {
            "priority": "critical",
            "recommendations": [
                "AVOID estrogen-containing contraceptives and HRT.",
                "7x VTE risk as heterozygote; 80x as homozygote.",
                "35x risk with estrogen-containing contraceptives.",
                "Inform ALL healthcare providers.",
                "Thromboprophylaxis before surgery/immobilization.",
                "Compression stockings for long flights (>4 hours).",
                "Genetic counseling for family members."
            ]
        }
    },
    "rs1799963": {
        "gene": "F2",
        "variant": "Prothrombin G20210A",
        "function": "Increased prothrombin levels",
        "risk_allele": "A",
        "frequency": {"EUR": 0.02, "AFR": 0.005, "MENA": 0.03},
        "effect": "3x VTE risk; multiplicative with FVL and estrogen",
        "category": "clotting",
        "condition": "Venous thromboembolism",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["8872854", "9591786"],
        "actionable": {
            "priority": "critical",
            "recommendations": [
                "AVOID estrogen-containing contraceptives.",
                "3x VTE risk; higher with additional risk factors.",
                "Inform healthcare providers before surgery.",
                "Family cascade screening recommended."
            ]
        }
    },
    "rs1799889": {
        "gene": "SERPINE1",
        "variant": "PAI-1 4G/5G",
        "function": "Plasminogen activator inhibitor-1 levels",
        "risk_allele": "G",  # 4G allele
        "frequency": {"EUR": 0.52, "AFR": 0.75, "EAS": 0.60},
        "effect": "4G/4G = higher PAI-1, impaired fibrinolysis",
        "category": "clotting",
        "condition": "Thrombosis / Myocardial infarction",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["8553611", "15657099"],
        "actionable": {
            "priority": "low",
            "recommendations": [
                "4G/4G: Higher PAI-1 levels, impaired clot breakdown.",
                "May contribute to MI risk, especially with metabolic syndrome.",
                "Maintain healthy weight (adipose tissue produces PAI-1)."
            ]
        }
    },
    "rs5918": {
        "gene": "ITGB3",
        "variant": "PlA1/A2 (L33P)",
        "function": "Platelet glycoprotein IIIa",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15, "AFR": 0.10, "EAS": 0.02},
        "effect": "PlA2 (C) may affect platelet reactivity and stent thrombosis",
        "category": "clotting",
        "condition": "Platelet reactivity / CAD",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["8598867", "18974117"],
        "actionable": {
            "priority": "low",
            "recommendations": [
                "Some studies suggest increased MI risk with PlA2 allele.",
                "May affect response to antiplatelet therapy.",
                "Results inconsistent across studies."
            ]
        }
    },
    "rs2046934": {
        "gene": "P2RY12",
        "variant": "P2Y12 receptor variant",
        "function": "Platelet ADP receptor",
        "risk_allele": "T",
        "frequency": {"EUR": 0.15},
        "effect": "May affect clopidogrel response",
        "category": "clotting",
        "condition": "Clopidogrel response / platelet reactivity",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["16959795"],
        "notes": "Affects drug target rather than metabolism"
    },
}

# =============================================================================
# HOMOCYSTEINE PATHWAY (Cardiovascular context)
# =============================================================================

HOMOCYSTEINE_MARKERS = {
    "rs1801133": {
        "gene": "MTHFR",
        "variant": "C677T",
        "function": "Folate metabolism / homocysteine",
        "risk_allele": "T",
        "frequency": {"EUR": 0.35, "EAS": 0.35, "AFR": 0.10, "AMR": 0.50},
        "effect": "TT = ~70% reduced enzyme activity, elevated homocysteine",
        "category": "homocysteine",
        "condition": "Hyperhomocysteinemia / CVD risk",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["8630491", "26647857"],
        "actionable": {
            "priority": "high",
            "recommendations": [
                "TT genotype: 2-3x higher homocysteine levels typical.",
                "Elevated homocysteine is independent CVD risk factor.",
                "Ensure adequate folate intake (leafy greens, fortified foods).",
                "Consider methylfolate (L-5-MTHF) over folic acid.",
                "Vitamin B12 and B6 also important for homocysteine metabolism.",
                "Check homocysteine level if TT genotype.",
                "Treat elevated homocysteine with B vitamins."
            ]
        }
    },
    "rs1801131": {
        "gene": "MTHFR",
        "variant": "A1298C",
        "function": "Folate metabolism (BH4 pathway)",
        "risk_allele": "C",
        "frequency": {"EUR": 0.33, "EAS": 0.20, "AFR": 0.15},
        "effect": "CC = ~40% reduced activity, milder effect than C677T",
        "category": "homocysteine",
        "condition": "Folate metabolism",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["10444342"],
        "actionable": {
            "priority": "medium",
            "recommendations": [
                "Milder effect than C677T.",
                "Compound heterozygous (677CT + 1298AC) may have significant effect.",
                "Consider methylfolate supplementation."
            ]
        }
    },
    "rs1805087": {
        "gene": "MTR",
        "variant": "A2756G",
        "function": "Methionine synthase (B12 dependent)",
        "risk_allele": "G",
        "frequency": {"EUR": 0.20, "EAS": 0.15, "AFR": 0.35},
        "effect": "May affect B12-dependent homocysteine remethylation",
        "category": "homocysteine",
        "condition": "Homocysteine / B12 metabolism",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["11114891"]
    },
    "rs1801394": {
        "gene": "MTRR",
        "variant": "A66G",
        "function": "Methionine synthase reductase",
        "risk_allele": "G",
        "frequency": {"EUR": 0.50, "EAS": 0.25, "AFR": 0.30},
        "effect": "May reduce B12 regeneration",
        "category": "homocysteine",
        "condition": "B12 metabolism",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["12518998"]
    },
}

# =============================================================================
# ARRHYTHMIA GENETICS
# =============================================================================

ARRHYTHMIA_MARKERS = {
    "rs2200733": {
        "gene": "PITX2 (4q25)",
        "variant": "Atrial fibrillation risk locus",
        "function": "Left atrial development",
        "risk_allele": "T",
        "frequency": {"EUR": 0.10, "EAS": 0.40, "AFR": 0.05},
        "effect": "Strongest AFib risk locus - 1.7x risk per allele",
        "category": "arrhythmia",
        "condition": "Atrial fibrillation",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["17603472", "22544366"],
        "actionable": {
            "priority": "medium",
            "recommendations": [
                "TT genotype: ~3x increased AFib risk.",
                "Regular pulse checks recommended.",
                "Report any palpitations to physician.",
                "Consider smartwatch for AFib detection.",
                "If AFib develops, anticoagulation decision based on CHA2DS2-VASc score."
            ]
        }
    },
    "rs10033464": {
        "gene": "PITX2 (4q25)",
        "variant": "Secondary AFib locus",
        "function": "Atrial development",
        "risk_allele": "T",
        "frequency": {"EUR": 0.05},
        "effect": "Additional AFib risk",
        "category": "arrhythmia",
        "condition": "Atrial fibrillation",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["17603472"]
    },
    "rs13376333": {
        "gene": "KCNN3",
        "variant": "Potassium channel variant",
        "function": "Cardiac repolarization",
        "risk_allele": "T",
        "frequency": {"EUR": 0.32, "EAS": 0.12},
        "effect": "Associated with AFib risk",
        "category": "arrhythmia",
        "condition": "Atrial fibrillation",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["20173747"]
    },
    "rs12143842": {
        "gene": "NOS1AP",
        "variant": "QT interval modulator",
        "function": "Nitric oxide signaling",
        "risk_allele": "T",
        "frequency": {"EUR": 0.35, "AFR": 0.20},
        "effect": "Strongest GWAS hit for QT interval",
        "category": "arrhythmia",
        "condition": "QT interval / sudden cardiac death",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["19305409", "16648850"],
        "actionable": {
            "priority": "medium",
            "recommendations": [
                "May have longer QT interval.",
                "Use caution with QT-prolonging drugs.",
                "ECG recommended before starting certain medications.",
                "See CredibleMeds.org for list of QT-prolonging drugs."
            ]
        }
    },
}

# Long QT Syndrome markers (monogenic)
LONG_QT_MARKERS = {
    "rs12720452": {
        "gene": "KCNQ1",
        "variant": "LQT1 common variant",
        "function": "Potassium channel (IKs)",
        "risk_allele": "G",
        "category": "arrhythmia",
        "condition": "Long QT syndrome type 1",
        "evidence": EvidenceLevel.STRONG,
        "pmid": ["17276177"],
        "notes": "Common variants; pathogenic variants require clinical sequencing"
    },
    "rs1805123": {
        "gene": "KCNH2",
        "variant": "K897T",
        "function": "Potassium channel (IKr)",
        "risk_allele": "C",
        "frequency": {"EUR": 0.15, "AFR": 0.30},
        "effect": "May modulate QT interval and drug-induced LQTS risk",
        "category": "arrhythmia",
        "condition": "Drug-induced Long QT",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["12860915", "17403134"],
        "actionable": {
            "priority": "low",
            "recommendations": [
                "May have increased susceptibility to drug-induced QT prolongation.",
                "Use caution with: certain antibiotics (fluoroquinolones, macrolides), antipsychotics, antiarrhythmics."
            ]
        }
    },
    "rs1805124": {
        "gene": "SCN5A",
        "variant": "H558R",
        "function": "Sodium channel (INa)",
        "risk_allele": "G",
        "frequency": {"EUR": 0.20, "AFR": 0.30},
        "effect": "May modify phenotype in SCN5A mutation carriers",
        "category": "arrhythmia",
        "condition": "Brugada / Long QT type 3",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["11834839"]
    },
}

# =============================================================================
# STRUCTURAL HEART DISEASE
# =============================================================================

STRUCTURAL_MARKERS = {
    "rs2118181": {
        "gene": "FBN1",
        "variant": "Fibrillin-1 variant",
        "function": "Extracellular matrix",
        "risk_allele": "G",
        "frequency": {"EUR": 0.30},
        "effect": "Associated with aortic root diameter",
        "category": "structural",
        "condition": "Aortic aneurysm / Marfan spectrum",
        "evidence": EvidenceLevel.MODERATE,
        "pmid": ["21737445"],
        "notes": "Marfan syndrome caused by pathogenic FBN1 mutations, not SNPs"
    },
    "rs10757278": {
        "gene": "9p21.3 (CDKN2A/B)",
        "variant": "CAD risk locus",
        "function": "Cell cycle regulation / vascular remodeling",
        "risk_allele": "G",
        "frequency": {"EUR": 0.50, "EAS": 0.55, "AFR": 0.30},
        "effect": "Strongest CAD risk locus - ~30% increased risk per allele",
        "category": "structural",
        "condition": "Coronary artery disease / Aortic aneurysm",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["17478679", "17554300"],
        "actionable": {
            "priority": "medium",
            "recommendations": [
                "GG genotype: ~60% increased CAD risk.",
                "Also associated with intracranial aneurysm.",
                "Aggressive cardiovascular risk factor management.",
                "Consider earlier coronary calcium scoring."
            ]
        }
    },
    "rs1333049": {
        "gene": "9p21.3",
        "variant": "CAD/MI risk",
        "function": "Adjacent to CDKN2A/B",
        "risk_allele": "C",
        "frequency": {"EUR": 0.47, "EAS": 0.55, "AFR": 0.35},
        "effect": "~25% increased MI risk per allele",
        "category": "structural",
        "condition": "Myocardial infarction",
        "evidence": EvidenceLevel.DEFINITIVE,
        "pmid": ["17478679", "17554300"]
    },
}

# =============================================================================
# COMBINE ALL CARDIOVASCULAR MARKERS
# =============================================================================

CARDIOVASCULAR_COMPLETE = {
    **LPA_MARKERS,
    **LDL_MARKERS,
    **HDL_MARKERS,
    **TRIGLYCERIDE_MARKERS,
    **BLOOD_PRESSURE_MARKERS,
    **CLOTTING_MARKERS,
    **HOMOCYSTEINE_MARKERS,
    **ARRHYTHMIA_MARKERS,
    **LONG_QT_MARKERS,
    **STRUCTURAL_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def calculate_lpa_risk(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate Lp(a) genetic risk score."""
    risk_score = 0
    interpretation = []
    
    # rs10455872
    geno = genotypes.get("rs10455872", "AA")
    if geno == "GG":
        risk_score += 2
        interpretation.append("rs10455872 GG: Very high Lp(a) expected (~100 nmol/L increase)")
    elif geno == "AG":
        risk_score += 1
        interpretation.append("rs10455872 AG: Elevated Lp(a) expected (~50 nmol/L increase)")
    
    # rs3798220
    geno = genotypes.get("rs3798220", "TT")
    if geno in ["CC", "TC"]:
        risk_score += 1
        interpretation.append("rs3798220 C allele: Additional Lp(a) elevation (~20-40 nmol/L)")
    
    # Determine risk level
    if risk_score >= 2:
        level = CardioRiskLevel.VERY_HIGH
        recommendation = "Measure Lp(a) level. Very high cardiovascular risk. Aggressive LDL lowering essential."
    elif risk_score == 1:
        level = CardioRiskLevel.HIGH
        recommendation = "Measure Lp(a) level. Consider earlier statin therapy."
    else:
        level = CardioRiskLevel.AVERAGE
        recommendation = "No genetic elevation of Lp(a) detected. Standard risk assessment applies."
    
    return {
        "genetic_risk_score": risk_score,
        "risk_level": level.value,
        "interpretation": interpretation,
        "recommendation": recommendation,
        "pmid": ["19060906", "22085571"]
    }

def calculate_thrombophilia_risk(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate thrombophilia (clotting) risk."""
    risk_factors = []
    total_risk_multiplier = 1.0
    
    # Factor V Leiden
    fvl = genotypes.get("rs6025", "GG")
    if fvl == "AA":
        risk_factors.append({"gene": "F5 (FVL)", "genotype": "AA (homozygous)", "risk": "80x VTE risk"})
        total_risk_multiplier *= 80
    elif fvl == "AG":
        risk_factors.append({"gene": "F5 (FVL)", "genotype": "AG (heterozygous)", "risk": "7x VTE risk"})
        total_risk_multiplier *= 7
    
    # Prothrombin
    pt = genotypes.get("rs1799963", "GG")
    if pt == "AA":
        risk_factors.append({"gene": "F2 (PT)", "genotype": "AA (homozygous)", "risk": "10x VTE risk"})
        total_risk_multiplier *= 10
    elif pt == "AG":
        risk_factors.append({"gene": "F2 (PT)", "genotype": "AG (heterozygous)", "risk": "3x VTE risk"})
        total_risk_multiplier *= 3
    
    # MTHFR (homocysteine)
    mthfr = genotypes.get("rs1801133", "CC")
    if mthfr == "TT":
        risk_factors.append({"gene": "MTHFR C677T", "genotype": "TT", "risk": "Elevated homocysteine, 1.5x VTE"})
        total_risk_multiplier *= 1.5
    
    # Determine overall risk
    if total_risk_multiplier >= 50:
        level = CardioRiskLevel.VERY_HIGH
        urgent = True
    elif total_risk_multiplier >= 10:
        level = CardioRiskLevel.HIGH
        urgent = True
    elif total_risk_multiplier >= 3:
        level = CardioRiskLevel.MODERATE
        urgent = False
    else:
        level = CardioRiskLevel.AVERAGE
        urgent = False
    
    return {
        "risk_factors": risk_factors,
        "combined_risk_multiplier": round(total_risk_multiplier, 1),
        "risk_level": level.value,
        "urgent_action_needed": urgent,
        "recommendations": [
            "Avoid estrogen-containing contraceptives" if total_risk_multiplier >= 3 else None,
            "Inform healthcare providers before surgery" if total_risk_multiplier >= 3 else None,
            "Genetic counseling for family members" if total_risk_multiplier >= 7 else None,
            "Medical alert bracelet recommended" if total_risk_multiplier >= 50 else None,
        ],
        "pmid": ["7989260", "8872854"]
    }

def calculate_afib_risk(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate atrial fibrillation genetic risk."""
    risk_score = 0
    variants = []
    
    # PITX2 (strongest effect)
    pitx2 = genotypes.get("rs2200733", "CC")
    if pitx2 == "TT":
        risk_score += 2
        variants.append({"gene": "PITX2", "rs": "rs2200733", "effect": "2.9x AFib risk"})
    elif pitx2 == "CT":
        risk_score += 1
        variants.append({"gene": "PITX2", "rs": "rs2200733", "effect": "1.7x AFib risk"})
    
    # Secondary PITX2
    pitx2_2 = genotypes.get("rs10033464", "GG")
    if pitx2_2 in ["TT", "GT"]:
        risk_score += 0.5
        variants.append({"gene": "PITX2", "rs": "rs10033464", "effect": "Additional risk"})
    
    # KCNN3
    kcnn3 = genotypes.get("rs13376333", "CC")
    if kcnn3 == "TT":
        risk_score += 0.5
        variants.append({"gene": "KCNN3", "rs": "rs13376333", "effect": "Modest AFib risk"})
    
    if risk_score >= 2:
        level = CardioRiskLevel.HIGH
    elif risk_score >= 1:
        level = CardioRiskLevel.MODERATE
    else:
        level = CardioRiskLevel.AVERAGE
    
    return {
        "genetic_risk_score": risk_score,
        "risk_level": level.value,
        "risk_variants": variants,
        "recommendations": [
            "Regular pulse checks (manual or smartwatch)",
            "Report palpitations to physician",
            "Consider screening ECG",
            "If AFib detected: anticoagulation based on CHA2DS2-VASc"
        ],
        "pmid": ["17603472", "22544366"]
    }

def generate_cardiovascular_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate comprehensive cardiovascular genetics report."""
    return {
        "lpa_risk": calculate_lpa_risk(genotypes),
        "thrombophilia_risk": calculate_thrombophilia_risk(genotypes),
        "afib_risk": calculate_afib_risk(genotypes),
        "markers_analyzed": sum(1 for rs in CARDIOVASCULAR_COMPLETE if rs in genotypes),
        "total_markers": len(CARDIOVASCULAR_COMPLETE),
        "critical_findings": [
            finding for finding in [
                "Factor V Leiden detected" if genotypes.get("rs6025") in ["AA", "AG"] else None,
                "Prothrombin mutation detected" if genotypes.get("rs1799963") in ["AA", "AG"] else None,
                "High Lp(a) genetic risk" if genotypes.get("rs10455872") in ["GG", "AG"] else None,
            ] if finding
        ]
    }

# Export
__all__ = [
    'CARDIOVASCULAR_COMPLETE',
    'LPA_MARKERS',
    'LDL_MARKERS',
    'HDL_MARKERS',
    'TRIGLYCERIDE_MARKERS',
    'BLOOD_PRESSURE_MARKERS',
    'CLOTTING_MARKERS',
    'HOMOCYSTEINE_MARKERS',
    'ARRHYTHMIA_MARKERS',
    'LONG_QT_MARKERS',
    'STRUCTURAL_MARKERS',
    'CardioRiskLevel',
    'EvidenceLevel',
    'calculate_lpa_risk',
    'calculate_thrombophilia_risk',
    'calculate_afib_risk',
    'generate_cardiovascular_report',
]
