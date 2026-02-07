"""
Personal Genomics v5.0.0

Comprehensive local DNA analysis with reference dataset integration.
Privacy-first: all data downloaded and processed locally.
"""

__version__ = "5.0.0"
__author__ = "OpenClaw"

from .datasets import (
    ThousandGenomes,
    GnomAD,
    ClinVar,
    PharmGKB,
    PGSCatalog,
    GWASCatalog,
    HGDP,
    SGDP,
    AncientDNA,
)

__all__ = [
    "ThousandGenomes",
    "GnomAD", 
    "ClinVar",
    "PharmGKB",
    "PGSCatalog",
    "GWASCatalog",
    "HGDP",
    "SGDP",
    "AncientDNA",
]
