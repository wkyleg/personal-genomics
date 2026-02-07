"""
Personal Genomics v5.0 Integration Module

This module integrates all v5.0 comprehensive marker expansions:
- pharmacogenomics_complete.py
- cardiovascular_complete.py
- supplement_protocol.py
- daily_optimization.py
- athletic_comprehensive.py
- nutrition_comprehensive.py
- mental_health_comprehensive.py
- skin_appearance.py
- longevity_comprehensive.py
- quirky_traits.py
- medical_special.py

Import this module to access all v5.0 features.
"""

from typing import Dict, List, Any, Optional
from pathlib import Path

# =============================================================================
# IMPORT ALL V5.0 MODULES
# =============================================================================

from .pharmacogenomics_complete import (
    PHARMACOGENOMICS_COMPLETE,
    CYP2D6_MARKERS, CYP2C19_MARKERS, CYP2C9_MARKERS, CYP3A_MARKERS,
    CYP1A2_MARKERS, SLCO1B1_MARKERS, VKORC1_MARKERS,
    DPYD_MARKERS, TPMT_MARKERS, NUDT15_MARKERS, UGT1A1_MARKERS,
    HLA_MARKERS, OPRM1_MARKERS, COMT_DRUG_MARKERS, ABCB1_MARKERS,
    NAT2_MARKERS, G6PD_MARKERS, ANESTHESIA_MARKERS, BCHE_MARKERS,
    MetabolizerStatus, ClinicalActionLevel,
    calculate_cyp2d6_activity_score, calculate_cyp2c19_status,
    get_critical_alerts, generate_pharmacogenomics_report
)

from .cardiovascular_complete import (
    CARDIOVASCULAR_COMPLETE,
    LPA_MARKERS, LDL_MARKERS, HDL_MARKERS, TRIGLYCERIDE_MARKERS,
    BLOOD_PRESSURE_MARKERS, CLOTTING_MARKERS, HOMOCYSTEINE_MARKERS,
    ARRHYTHMIA_MARKERS, LONG_QT_MARKERS, STRUCTURAL_MARKERS,
    CardioRiskLevel, EvidenceLevel,
    calculate_lpa_risk, calculate_thrombophilia_risk,
    calculate_afib_risk, generate_cardiovascular_report
)

from .supplement_protocol import (
    SUPPLEMENT_MARKERS,
    METHYLATION_MARKERS, VITAMIN_D_MARKERS, VITAMIN_A_MARKERS,
    VITAMIN_B12_MARKERS, VITAMIN_E_MARKERS, OMEGA3_MARKERS,
    IRON_MARKERS, MAGNESIUM_MARKERS, ZINC_MARKERS,
    ANTIOXIDANT_MARKERS, COQ10_MARKERS, CHOLINE_MARKERS,
    SupplementPriority,
    generate_supplement_protocol
)

from .daily_optimization import (
    DAILY_OPTIMIZATION_MARKERS,
    CHRONOTYPE_MARKERS, CAFFEINE_MARKERS,
    EXERCISE_TIMING_MARKERS, MEAL_TIMING_MARKERS,
    Chronotype, CaffeineMetabolism, ExerciseTimingType,
    calculate_chronotype_score, calculate_caffeine_profile,
    calculate_optimal_exercise_time, calculate_optimal_meal_timing,
    generate_daily_optimization_report
)

from .athletic_comprehensive import (
    ATHLETIC_COMPREHENSIVE_MARKERS,
    POWER_MARKERS, ENDURANCE_MARKERS, LACTATE_MARKERS,
    RECOVERY_MARKERS, TENDON_MARKERS, BONE_MARKERS,
    CREATINE_MARKERS, CAFFEINE_ERGOGENIC_MARKERS,
    HEAT_MARKERS, ALTITUDE_MARKERS, CONCUSSION_MARKERS, HYDRATION_MARKERS,
    AthleticProfile, RecoverySpeed, InjuryRiskLevel,
    calculate_athletic_profile, calculate_injury_risk,
    calculate_recovery_profile, generate_athletic_report
)

from .nutrition_comprehensive import (
    NUTRITION_COMPREHENSIVE_MARKERS,
    FAT_METABOLISM_MARKERS, APOE_MARKERS, APOE_DIET_RECOMMENDATIONS,
    CARB_MARKERS, SODIUM_MARKERS, HISTAMINE_MARKERS,
    OXALATE_MARKERS, CELIAC_MARKERS, LACTOSE_MARKERS,
    FRUCTOSE_MARKERS, ALCOHOL_MARKERS, TASTE_MARKERS,
    ToleranceLevel, DietaryRecommendation,
    determine_apoe_genotype, generate_nutrition_report
)

from .mental_health_comprehensive import (
    MENTAL_HEALTH_COMPREHENSIVE_MARKERS, MENTAL_HEALTH_DISCLAIMER,
    DEPRESSION_MARKERS, ANTIDEPRESSANT_MARKERS, ANXIETY_MARKERS,
    ADHD_MARKERS, BIPOLAR_MARKERS, STRESS_MARKERS,
    SOCIAL_MARKERS, AGGRESSION_MARKERS, PERSONALITY_MARKERS, ADDICTION_MARKERS,
    RiskLevel, StressResponseType, EvidenceStrength,
    determine_stress_type, assess_depression_vulnerability,
    assess_addiction_vulnerability, generate_mental_health_report
)

from .skin_appearance import (
    SKIN_APPEARANCE_MARKERS,
    PIGMENTATION_MARKERS, MC1R_MARKERS, MELANOMA_MARKERS,
    SKIN_AGING_MARKERS, SKIN_CONDITION_MARKERS, HAIR_MARKERS,
    SkinType, MelanomaRisk,
    calculate_skin_type, calculate_melanoma_risk,
    predict_eye_color, predict_hair_loss_risk,
    generate_skin_appearance_report
)

from .longevity_comprehensive import (
    LONGEVITY_COMPREHENSIVE_MARKERS,
    FOXO3_MARKERS, APOE_LONGEVITY, TERT_MARKERS,
    CETP_MARKERS, IL6_LONGEVITY, KLOTHO_MARKERS, SIRT_MARKERS,
    LongevityPotential,
    calculate_longevity_score, generate_longevity_report
)

from .quirky_traits import (
    QUIRKY_TRAITS_MARKERS,
    ABCC11_MARKERS, ASPARAGUS_MARKERS, PHOTIC_SNEEZE_MARKERS,
    CILANTRO_MARKERS, MOSQUITO_MARKERS, MOTION_SICKNESS_MARKERS,
    ACROPHOBIA_MARKERS, OTHER_FUN_MARKERS,
    generate_quirky_traits_report
)

from .medical_special import (
    MEDICAL_SPECIAL_MARKERS,
    ANESTHESIA_MH_MARKERS, BCHE_MARKERS as BCHE_ANESTHESIA_MARKERS,
    ABO_MARKERS, RH_MARKERS, INFECTION_RESISTANCE_MARKERS,
    CELIAC_HLA_COMPLETE, VACCINE_RESPONSE_MARKERS,
    AnesthesiaRisk, BloodType,
    check_anesthesia_alerts, infer_blood_type,
    check_celiac_hla, check_infection_resistance,
    generate_medical_special_report
)

# =============================================================================
# COMBINED V5.0 MARKER DATABASE
# =============================================================================

V5_ALL_MARKERS = {
    **PHARMACOGENOMICS_COMPLETE,
    **CARDIOVASCULAR_COMPLETE,
    **SUPPLEMENT_MARKERS,
    **DAILY_OPTIMIZATION_MARKERS,
    **ATHLETIC_COMPREHENSIVE_MARKERS,
    **NUTRITION_COMPREHENSIVE_MARKERS,
    **MENTAL_HEALTH_COMPREHENSIVE_MARKERS,
    **SKIN_APPEARANCE_MARKERS,
    **LONGEVITY_COMPREHENSIVE_MARKERS,
    **QUIRKY_TRAITS_MARKERS,
    **MEDICAL_SPECIAL_MARKERS,
}

# =============================================================================
# MASTER ANALYSIS FUNCTION
# =============================================================================

def parse_dna_file(filepath: str) -> Dict[str, str]:
    """Parse DNA file (23andMe or AncestryDNA format) into genotypes dict."""
    genotypes = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) >= 4:
                rsid = parts[0]
                # Combine alleles
                allele1 = parts[3] if len(parts) > 3 else ''
                allele2 = parts[4] if len(parts) > 4 else ''
                
                if allele1 and allele2:
                    genotypes[rsid] = allele1 + allele2
                elif allele1:
                    genotypes[rsid] = allele1
    
    return genotypes

def generate_comprehensive_v5_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate complete v5.0 genomics report."""
    
    # Count markers found
    markers_found = sum(1 for rs in V5_ALL_MARKERS if rs in genotypes)
    total_markers = len(V5_ALL_MARKERS)
    
    report = {
        "version": "5.0",
        "markers_analyzed": markers_found,
        "total_v5_markers": total_markers,
        "coverage_percent": round(100 * markers_found / total_markers, 1) if total_markers > 0 else 0,
        
        # Critical medical alerts first
        "critical_alerts": {
            "pharmacogenomics": get_critical_alerts(genotypes),
            "anesthesia": check_anesthesia_alerts(genotypes),
        },
        
        # Core medical categories
        "pharmacogenomics": generate_pharmacogenomics_report(genotypes),
        "cardiovascular": generate_cardiovascular_report(genotypes),
        "supplement_protocol": generate_supplement_protocol(genotypes),
        "medical_special": generate_medical_special_report(genotypes),
        
        # Lifestyle optimization
        "daily_optimization": generate_daily_optimization_report(genotypes),
        "athletic": generate_athletic_report(genotypes),
        "nutrition": generate_nutrition_report(genotypes),
        
        # Health and traits
        "mental_health": generate_mental_health_report(genotypes),
        "skin_appearance": generate_skin_appearance_report(genotypes),
        "longevity": generate_longevity_report(genotypes),
        
        # Fun stuff
        "quirky_traits": generate_quirky_traits_report(genotypes),
    }
    
    # Summary of critical findings
    all_critical = []
    
    if report["critical_alerts"]["pharmacogenomics"]:
        all_critical.extend([
            f"⚠️ {a['gene']}: {a['drug']} - {a['message']}"
            for a in report["critical_alerts"]["pharmacogenomics"]
        ])
    
    if report["critical_alerts"]["anesthesia"]["has_critical"]:
        all_critical.extend([
            f"🚨 {a['gene']}: {a['condition']} - {a['message']}"
            for a in report["critical_alerts"]["anesthesia"]["alerts"]
            if a["severity"] == "CRITICAL"
        ])
    
    report["critical_findings_summary"] = all_critical
    
    return report

def get_v5_marker_counts() -> Dict[str, int]:
    """Get marker counts by category for v5.0."""
    return {
        "pharmacogenomics": len(PHARMACOGENOMICS_COMPLETE),
        "cardiovascular": len(CARDIOVASCULAR_COMPLETE),
        "supplement": len(SUPPLEMENT_MARKERS),
        "daily_optimization": len(DAILY_OPTIMIZATION_MARKERS),
        "athletic": len(ATHLETIC_COMPREHENSIVE_MARKERS),
        "nutrition": len(NUTRITION_COMPREHENSIVE_MARKERS),
        "mental_health": len(MENTAL_HEALTH_COMPREHENSIVE_MARKERS),
        "skin_appearance": len(SKIN_APPEARANCE_MARKERS),
        "longevity": len(LONGEVITY_COMPREHENSIVE_MARKERS),
        "quirky_traits": len(QUIRKY_TRAITS_MARKERS),
        "medical_special": len(MEDICAL_SPECIAL_MARKERS),
        "total_unique": len(V5_ALL_MARKERS),
    }

# =============================================================================
# EXPORTS
# =============================================================================

__all__ = [
    # Master collections
    'V5_ALL_MARKERS',
    'get_v5_marker_counts',
    
    # Master analysis
    'parse_dna_file',
    'generate_comprehensive_v5_report',
    
    # Pharmacogenomics
    'PHARMACOGENOMICS_COMPLETE',
    'get_critical_alerts',
    'generate_pharmacogenomics_report',
    'MetabolizerStatus',
    'ClinicalActionLevel',
    
    # Cardiovascular
    'CARDIOVASCULAR_COMPLETE',
    'generate_cardiovascular_report',
    'CardioRiskLevel',
    
    # Supplements
    'SUPPLEMENT_MARKERS',
    'generate_supplement_protocol',
    
    # Daily optimization
    'DAILY_OPTIMIZATION_MARKERS',
    'generate_daily_optimization_report',
    'Chronotype',
    'CaffeineMetabolism',
    
    # Athletic
    'ATHLETIC_COMPREHENSIVE_MARKERS',
    'generate_athletic_report',
    'AthleticProfile',
    
    # Nutrition
    'NUTRITION_COMPREHENSIVE_MARKERS',
    'generate_nutrition_report',
    'determine_apoe_genotype',
    
    # Mental Health
    'MENTAL_HEALTH_COMPREHENSIVE_MARKERS',
    'generate_mental_health_report',
    'StressResponseType',
    
    # Skin/Appearance
    'SKIN_APPEARANCE_MARKERS',
    'generate_skin_appearance_report',
    
    # Longevity
    'LONGEVITY_COMPREHENSIVE_MARKERS',
    'generate_longevity_report',
    
    # Quirky traits
    'QUIRKY_TRAITS_MARKERS',
    'generate_quirky_traits_report',
    
    # Medical special
    'MEDICAL_SPECIAL_MARKERS',
    'generate_medical_special_report',
    'check_anesthesia_alerts',
    'check_celiac_hla',
    'check_infection_resistance',
]
