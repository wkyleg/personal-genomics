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
from .polygenic_scores import PRS_WEIGHTS, PRS_CONDITIONS
from .carrier_status import CARRIER_MARKERS
from .health_risks import HEALTH_RISK_MARKERS
from .traits import TRAIT_MARKERS
from .nutrition import NUTRITION_MARKERS
from .fitness import FITNESS_MARKERS
from .neurogenetics import NEURO_MARKERS
from .longevity import LONGEVITY_MARKERS
from .immunity import IMMUNITY_MARKERS

__all__ = [
    'PHARMACOGENOMICS_MARKERS',
    'DRUG_INTERACTIONS', 
    'PRS_WEIGHTS',
    'PRS_CONDITIONS',
    'CARRIER_MARKERS',
    'HEALTH_RISK_MARKERS',
    'TRAIT_MARKERS',
    'NUTRITION_MARKERS',
    'FITNESS_MARKERS',
    'NEURO_MARKERS',
    'LONGEVITY_MARKERS',
    'IMMUNITY_MARKERS'
]
