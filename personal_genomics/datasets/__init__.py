"""
Genomics Reference Datasets

This package provides access to major population genetics and clinical databases:

POPULATION/ANCESTRY:
- 1000 Genomes Project (26 populations, ~2,500 individuals)
- HGDP - Human Genome Diversity Project (51 populations, ~1,000 individuals)
- SGDP - Simons Genome Diversity Project (142 populations, 300 individuals)
- Ancient DNA (WHG, Neolithic, Yamnaya, Neanderthal markers)

HEALTH/CLINICAL:
- gnomAD - Genome Aggregation Database (allele frequencies)
- ClinVar - Clinical variant interpretations
- PharmGKB - Pharmacogenomics database
- PGS Catalog - Polygenic risk scores
- GWAS Catalog - Genome-wide associations

All datasets are downloaded and cached locally for privacy-first analysis.
No data is sent to external servers during analysis.
"""

from .base import (
    BaseDataset,
    SQLiteDataset,
    DatasetVersion,
    VariantInfo,
    PopulationFrequency,
    DATASETS_BASE_PATH,
    THOUSAND_GENOMES_POPULATIONS,
    SUPERPOPULATIONS,
    get_all_datasets,
    download_all_datasets,
    get_dataset_status,
)

from .thousand_genomes import (
    ThousandGenomes,
    get_1kg_frequencies,
    compare_to_1kg_populations,
)

from .hgdp import (
    HGDP,
    HGDP_POPULATIONS,
    HGDP_REGIONS,
    compare_to_hgdp,
)

from .sgdp import (
    SGDPDataset,
    SGDPPopulation,
    SGDP_POPULATIONS,
    ALL_SGDP_POPULATIONS,
    get_sgdp_info,
    list_european_populations,
    list_middle_eastern_populations,
)
# Alias for compatibility
SGDP = SGDPDataset

from .ancient_dna import (
    AncientDNADataset,
    AncientMarker,
    AncestralSignal,
    ANCIENT_MARKERS,
    MARKERS_BY_POPULATION,
    get_ancient_markers,
    analyze_ancient_ancestry,
    get_european_ancestry_narrative,
)
# Alias for compatibility
AncientDNA = AncientDNADataset

from .gnomad import (
    GnomAD,
    get_gnomad_frequency,
    is_rare_in_gnomad,
)

from .clinvar import (
    ClinVar,
    ClinVarVariant,
    get_clinical_significance,
    is_pathogenic,
)

from .pharmgkb import (
    PharmGKB,
    DrugGeneInteraction,
    DosingGuideline,
    get_drug_recommendations,
    check_medication_safety,
)

from .pgs_catalog import (
    PGSCatalog,
    PolygenticScore,
    PRSResult,
    calculate_disease_risk,
)

from .gwas_catalog import (
    GWASCatalog,
    GWASAssociation,
    get_trait_associations,
)

__all__ = [
    # Base classes
    "BaseDataset",
    "SQLiteDataset", 
    "DatasetVersion",
    "VariantInfo",
    "PopulationFrequency",
    "DATASETS_BASE_PATH",
    "get_all_datasets",
    "download_all_datasets",
    "get_dataset_status",
    
    # 1000 Genomes
    "ThousandGenomes",
    "THOUSAND_GENOMES_POPULATIONS",
    "SUPERPOPULATIONS",
    "get_1kg_frequencies",
    "compare_to_1kg_populations",
    
    # HGDP
    "HGDP",
    "HGDP_POPULATIONS",
    "HGDP_REGIONS",
    "compare_to_hgdp",
    
    # SGDP
    "SGDPDataset",
    "SGDP",  # Alias
    "SGDPPopulation",
    "SGDP_POPULATIONS",
    "ALL_SGDP_POPULATIONS",
    "get_sgdp_info",
    "list_european_populations",
    "list_middle_eastern_populations",
    
    # Ancient DNA
    "AncientDNADataset",
    "AncientDNA",  # Alias
    "AncientMarker",
    "AncestralSignal",
    "ANCIENT_MARKERS",
    "MARKERS_BY_POPULATION",
    "get_ancient_markers",
    "analyze_ancient_ancestry",
    "get_european_ancestry_narrative",
    
    # gnomAD
    "GnomAD",
    "get_gnomad_frequency",
    "is_rare_in_gnomad",
    
    # ClinVar
    "ClinVar",
    "ClinVarVariant",
    "get_clinical_significance",
    "is_pathogenic",
    
    # PharmGKB
    "PharmGKB",
    "DrugGeneInteraction",
    "DosingGuideline",
    "get_drug_recommendations",
    "check_medication_safety",
    
    # PGS Catalog
    "PGSCatalog",
    "PolygenticScore",
    "PRSResult",
    "calculate_disease_risk",
    
    # GWAS Catalog
    "GWASCatalog",
    "GWASAssociation",
    "get_trait_associations",
]
