"""
Personal Genomics Marker Database v4.0
Comprehensive collection of validated genetic markers from:
- PharmGKB (pharmacogenomics)
- ClinVar (clinical variants)
- NHGRI-EBI GWAS Catalog
- dbSNP
- CPIC (Clinical Pharmacogenetics Implementation Consortium)
- OMIM (rare diseases)
- PhyloTree / ISOGG (haplogroups)
- 1000 Genomes (ancestry)

All markers are validated, published, and appropriate for consumer genomics.

Categories:
1. Pharmacogenomics - Drug metabolism and response
2. Polygenic Risk Scores - Disease risk prediction
3. Carrier Status - Recessive disease carriers
4. Health Risks - Disease susceptibility
5. Traits - Physical and behavioral traits
6. Nutrition - Nutrigenomics
7. Fitness - Athletic performance and injury
8. Neurogenetics - Cognition and behavior
9. Longevity - Aging markers
10. Immunity - HLA and immune function
11. Ancestry - Population informative markers
12. Rare Diseases - Rare genetic conditions
13. Mental Health - Psychiatric genetics
14. Dermatology - Skin and hair
15. Vision & Hearing - Sensory genetics
16. Fertility - Reproductive health
17. Haplogroups - Maternal/paternal lineage (NEW v4.0)
18. Ancestry Composition - Population admixture (NEW v4.0)
19. Cancer Panel - Hereditary cancer markers (NEW v4.0)
20. Autoimmune HLA - HLA disease associations (NEW v4.0)
21. Pain Sensitivity - Pain perception/opioid response (NEW v4.0)
"""

from .pharmacogenomics import PHARMACOGENOMICS_MARKERS, DRUG_INTERACTIONS
from .pharmacogenomics_extended import PHARMACOGENOMICS_EXTENDED
from .polygenic_scores import PRS_WEIGHTS, PRS_CONDITIONS, calculate_prs
from .prs_extended import PRS_EXTENDED
from .carrier_status import CARRIER_MARKERS, CARRIER_SCREENING_PANELS
from .carrier_extended import CARRIER_EXTENDED
from .health_risks import HEALTH_RISK_MARKERS
from .health_extended import HEALTH_EXTENDED
from .traits import TRAIT_MARKERS
from .traits_extended import TRAITS_EXTENDED
from .nutrition import NUTRITION_MARKERS
from .fitness import FITNESS_MARKERS
from .neurogenetics import NEURO_MARKERS
from .longevity import LONGEVITY_MARKERS
from .immunity import IMMUNITY_MARKERS, HLA_DRUG_ALERTS
from .ancestry import ANCESTRY_MARKERS, POPULATION_CODES

# v3.0 categories
from .rare_diseases import RARE_DISEASE_MARKERS, RARE_DISEASE_PANELS
from .mental_health import MENTAL_HEALTH_MARKERS, MENTAL_HEALTH_NOTES
from .dermatology import DERMATOLOGY_MARKERS, DERMATOLOGY_SUMMARY
from .vision_hearing import VISION_MARKERS, HEARING_MARKERS, VISION_HEARING_MARKERS, CLINICAL_NOTES as SENSORY_NOTES
from .fertility import FERTILITY_MARKERS, REPRODUCTIVE_NOTES

# NEW v4.0 categories
from .haplogroups import (
    MTDNA_MARKERS, YCHROMOSOME_MARKERS, HAPLOGROUP_HISTORY,
    determine_mtdna_haplogroup, determine_y_haplogroup, analyze_haplogroups
)
from .ancestry_composition import (
    ANCESTRY_INFORMATIVE_MARKERS, POPULATION_DESCRIPTIONS,
    estimate_ancestry, detect_admixture, get_ancestry_summary
)
from .cancer_panel import (
    BRCA1_MARKERS, BRCA2_MARKERS, LYNCH_SYNDROME_MARKERS, OTHER_CANCER_MARKERS,
    HEREDITARY_CANCER_MARKERS, CANCER_SCREENING_PANELS, analyze_cancer_panel
)
from .autoimmune_hla import (
    CELIAC_HLA_MARKERS, TYPE1_DIABETES_HLA_MARKERS, ANKYLOSING_SPONDYLITIS_MARKERS,
    RHEUMATOID_ARTHRITIS_MARKERS, LUPUS_MARKERS, OTHER_AUTOIMMUNE_MARKERS,
    AUTOIMMUNE_HLA_MARKERS, AUTOIMMUNE_CONDITIONS, analyze_autoimmune_risk
)
from .pain_sensitivity import (
    COMT_MARKERS, OPRM1_MARKERS, SCN9A_MARKERS, TRPV1_MARKERS,
    OPIOID_METABOLISM_MARKERS, MIGRAINE_MARKERS, NSAID_RESPONSE_MARKERS,
    PAIN_SENSITIVITY_MARKERS, analyze_pain_sensitivity
)

# Merge extended markers into main dictionaries
def _merge_markers():
    """Merge extended markers into main dictionaries."""
    # Merge pharmacogenomics
    all_pharma = {**PHARMACOGENOMICS_MARKERS, **PHARMACOGENOMICS_EXTENDED}
    
    # Merge health risks (including autoimmune HLA)
    all_health = {**HEALTH_RISK_MARKERS, **HEALTH_EXTENDED, **AUTOIMMUNE_HLA_MARKERS}
    
    # Merge traits (including dermatology traits)
    all_traits = {**TRAIT_MARKERS, **TRAITS_EXTENDED}
    
    # Merge PRS
    all_prs = {**PRS_WEIGHTS, **PRS_EXTENDED}
    
    # Merge carrier (including rare diseases and cancer panel)
    all_carrier = {**CARRIER_MARKERS, **CARRIER_EXTENDED, **RARE_DISEASE_MARKERS, **HEREDITARY_CANCER_MARKERS}
    
    # Merge ancestry (including AIMs and haplogroups)
    all_ancestry = {**ANCESTRY_MARKERS, **ANCESTRY_INFORMATIVE_MARKERS, **MTDNA_MARKERS, **YCHROMOSOME_MARKERS}
    
    return all_pharma, all_health, all_traits, all_prs, all_carrier, all_ancestry

ALL_PHARMACOGENOMICS, ALL_HEALTH_RISKS, ALL_TRAITS, ALL_PRS, ALL_CARRIER, ALL_ANCESTRY = _merge_markers()

# Combined marker count
def get_marker_counts():
    """Get counts of markers in each category."""
    return {
        "pharmacogenomics": len(ALL_PHARMACOGENOMICS),
        "prs_weights": len(ALL_PRS),
        "carrier_status": len(ALL_CARRIER),
        "health_risks": len(ALL_HEALTH_RISKS),
        "traits": len(ALL_TRAITS),
        "nutrition": len(NUTRITION_MARKERS),
        "fitness": len(FITNESS_MARKERS),
        "neurogenetics": len(NEURO_MARKERS),
        "longevity": len(LONGEVITY_MARKERS),
        "immunity": len(IMMUNITY_MARKERS),
        "ancestry": len(ALL_ANCESTRY),
        "rare_diseases": len(RARE_DISEASE_MARKERS),
        "mental_health": len(MENTAL_HEALTH_MARKERS),
        "dermatology": len(DERMATOLOGY_MARKERS),
        "vision_hearing": len(VISION_HEARING_MARKERS),
        "fertility": len(FERTILITY_MARKERS),
        # New v4.0 categories
        "haplogroups": len(MTDNA_MARKERS) + len(YCHROMOSOME_MARKERS),
        "ancestry_aims": len(ANCESTRY_INFORMATIVE_MARKERS),
        "hereditary_cancer": len(HEREDITARY_CANCER_MARKERS),
        "autoimmune_hla": len(AUTOIMMUNE_HLA_MARKERS),
        "pain_sensitivity": len(PAIN_SENSITIVITY_MARKERS),
        "total": (len(ALL_PHARMACOGENOMICS) + len(ALL_PRS) + len(ALL_CARRIER) +
                 len(ALL_HEALTH_RISKS) + len(ALL_TRAITS) + len(NUTRITION_MARKERS) +
                 len(FITNESS_MARKERS) + len(NEURO_MARKERS) + len(LONGEVITY_MARKERS) +
                 len(IMMUNITY_MARKERS) + len(ALL_ANCESTRY) +
                 len(MENTAL_HEALTH_MARKERS) + len(DERMATOLOGY_MARKERS) +
                 len(VISION_HEARING_MARKERS) + len(FERTILITY_MARKERS) +
                 len(PAIN_SENSITIVITY_MARKERS))
    }

__all__ = [
    # Core pharmacogenomics
    'PHARMACOGENOMICS_MARKERS',
    'PHARMACOGENOMICS_EXTENDED',
    'ALL_PHARMACOGENOMICS',
    'DRUG_INTERACTIONS',
    
    # PRS
    'PRS_WEIGHTS',
    'PRS_EXTENDED',
    'ALL_PRS',
    'PRS_CONDITIONS',
    'calculate_prs',
    
    # Carrier status
    'CARRIER_MARKERS',
    'CARRIER_EXTENDED',
    'ALL_CARRIER',
    'CARRIER_SCREENING_PANELS',
    
    # Health risks
    'HEALTH_RISK_MARKERS',
    'HEALTH_EXTENDED',
    'ALL_HEALTH_RISKS',
    
    # Traits
    'TRAIT_MARKERS',
    'TRAITS_EXTENDED',
    'ALL_TRAITS',
    
    # Specialty categories
    'NUTRITION_MARKERS',
    'FITNESS_MARKERS',
    'NEURO_MARKERS',
    'LONGEVITY_MARKERS',
    'IMMUNITY_MARKERS',
    'HLA_DRUG_ALERTS',
    'ANCESTRY_MARKERS',
    'POPULATION_CODES',
    
    # v3.0 categories
    'RARE_DISEASE_MARKERS',
    'RARE_DISEASE_PANELS',
    'MENTAL_HEALTH_MARKERS',
    'MENTAL_HEALTH_NOTES',
    'DERMATOLOGY_MARKERS',
    'DERMATOLOGY_SUMMARY',
    'VISION_MARKERS',
    'HEARING_MARKERS',
    'VISION_HEARING_MARKERS',
    'SENSORY_NOTES',
    'FERTILITY_MARKERS',
    'REPRODUCTIVE_NOTES',
    
    # NEW v4.0 - Haplogroups
    'MTDNA_MARKERS',
    'YCHROMOSOME_MARKERS',
    'HAPLOGROUP_HISTORY',
    'determine_mtdna_haplogroup',
    'determine_y_haplogroup',
    'analyze_haplogroups',
    
    # NEW v4.0 - Ancestry Composition
    'ANCESTRY_INFORMATIVE_MARKERS',
    'POPULATION_DESCRIPTIONS',
    'ALL_ANCESTRY',
    'estimate_ancestry',
    'detect_admixture',
    'get_ancestry_summary',
    
    # NEW v4.0 - Hereditary Cancer
    'BRCA1_MARKERS',
    'BRCA2_MARKERS',
    'LYNCH_SYNDROME_MARKERS',
    'OTHER_CANCER_MARKERS',
    'HEREDITARY_CANCER_MARKERS',
    'CANCER_SCREENING_PANELS',
    'analyze_cancer_panel',
    
    # NEW v4.0 - Autoimmune HLA
    'CELIAC_HLA_MARKERS',
    'TYPE1_DIABETES_HLA_MARKERS',
    'ANKYLOSING_SPONDYLITIS_MARKERS',
    'RHEUMATOID_ARTHRITIS_MARKERS',
    'LUPUS_MARKERS',
    'OTHER_AUTOIMMUNE_MARKERS',
    'AUTOIMMUNE_HLA_MARKERS',
    'AUTOIMMUNE_CONDITIONS',
    'analyze_autoimmune_risk',
    
    # NEW v4.0 - Pain Sensitivity
    'COMT_MARKERS',
    'OPRM1_MARKERS',
    'SCN9A_MARKERS',
    'TRPV1_MARKERS',
    'OPIOID_METABOLISM_MARKERS',
    'MIGRAINE_MARKERS',
    'NSAID_RESPONSE_MARKERS',
    'PAIN_SENSITIVITY_MARKERS',
    'analyze_pain_sensitivity',
    
    # Utility
    'get_marker_counts'
]
