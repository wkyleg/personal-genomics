"""
Personal Genomics Marker Database
Comprehensive collection of validated genetic markers from:
- PharmGKB (pharmacogenomics)
- ClinVar (clinical variants)
- NHGRI-EBI GWAS Catalog
- dbSNP
- SNPedia

All markers are validated, published, and appropriate for consumer genomics.
"""

from .pharmacogenomics import PHARMACOGENOMICS_MARKERS, DRUG_INTERACTIONS
from .pharmacogenomics_extended import PHARMACOGENOMICS_EXTENDED
from .polygenic_scores import PRS_WEIGHTS, PRS_CONDITIONS
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

# Merge extended markers into main dictionaries
def _merge_markers():
    """Merge extended markers into main dictionaries."""
    # Merge pharmacogenomics
    all_pharma = {**PHARMACOGENOMICS_MARKERS, **PHARMACOGENOMICS_EXTENDED}
    
    # Merge health risks
    all_health = {**HEALTH_RISK_MARKERS, **HEALTH_EXTENDED}
    
    # Merge traits
    all_traits = {**TRAIT_MARKERS, **TRAITS_EXTENDED}
    
    # Merge PRS
    all_prs = {**PRS_WEIGHTS, **PRS_EXTENDED}
    
    # Merge carrier
    all_carrier = {**CARRIER_MARKERS, **CARRIER_EXTENDED}
    
    return all_pharma, all_health, all_traits, all_prs, all_carrier

ALL_PHARMACOGENOMICS, ALL_HEALTH_RISKS, ALL_TRAITS, ALL_PRS, ALL_CARRIER = _merge_markers()

# Combined marker count
def get_marker_counts():
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
        "ancestry": len(ANCESTRY_MARKERS),
        "total": (len(ALL_PHARMACOGENOMICS) + len(ALL_PRS) + len(ALL_CARRIER) +
                 len(ALL_HEALTH_RISKS) + len(ALL_TRAITS) + len(NUTRITION_MARKERS) +
                 len(FITNESS_MARKERS) + len(NEURO_MARKERS) + len(LONGEVITY_MARKERS) +
                 len(IMMUNITY_MARKERS) + len(ANCESTRY_MARKERS))
    }

__all__ = [
    'PHARMACOGENOMICS_MARKERS',
    'PHARMACOGENOMICS_EXTENDED',
    'ALL_PHARMACOGENOMICS',
    'DRUG_INTERACTIONS',
    'PRS_WEIGHTS',
    'PRS_EXTENDED',
    'ALL_PRS',
    'PRS_CONDITIONS',
    'CARRIER_MARKERS',
    'CARRIER_EXTENDED',
    'ALL_CARRIER',
    'CARRIER_SCREENING_PANELS',
    'HEALTH_RISK_MARKERS',
    'HEALTH_EXTENDED',
    'ALL_HEALTH_RISKS',
    'TRAIT_MARKERS',
    'TRAITS_EXTENDED',
    'ALL_TRAITS',
    'NUTRITION_MARKERS',
    'FITNESS_MARKERS',
    'NEURO_MARKERS',
    'LONGEVITY_MARKERS',
    'IMMUNITY_MARKERS',
    'HLA_DRUG_ALERTS',
    'ANCESTRY_MARKERS',
    'POPULATION_CODES',
    'get_marker_counts'
]
