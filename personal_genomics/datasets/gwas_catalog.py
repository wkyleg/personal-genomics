"""
GWAS Catalog (NHGRI-EBI) Integration

Provides genome-wide association study results for trait associations.
Links variants to traits with supporting literature.

Data source: https://www.ebi.ac.uk/gwas/
"""

from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
import json
import logging
import sqlite3

from .base import (
    SQLiteDataset,
    DatasetVersion,
    VariantInfo,
)

logger = logging.getLogger(__name__)


@dataclass
class GWASAssociation:
    """A GWAS association between a variant and a trait."""
    rsid: str
    trait: str
    trait_ontology: str
    p_value: float
    odds_ratio: Optional[float]
    beta: Optional[float]
    ci_lower: Optional[float]
    ci_upper: Optional[float]
    risk_allele: str
    gene: str
    pmid: str
    study_sample_size: int
    ancestry: str
    
    @property
    def effect_direction(self) -> str:
        """Get effect direction (risk/protective)."""
        if self.odds_ratio:
            return "risk" if self.odds_ratio > 1 else "protective"
        elif self.beta:
            return "increasing" if self.beta > 0 else "decreasing"
        return "unknown"


class GWASCatalog(SQLiteDataset):
    """
    GWAS Catalog integration.
    
    Provides:
    - Trait-variant associations
    - Effect sizes and p-values
    - Literature references
    
    Uses curated subset of high-confidence associations.
    """
    
    name = "gwas_catalog"
    version = "2024-02"
    description = "NHGRI-EBI GWAS Catalog"
    source_url = "https://www.ebi.ac.uk/gwas/"
    
    SCHEMA = """
    CREATE TABLE IF NOT EXISTS associations (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rsid TEXT,
        trait TEXT,
        trait_ontology TEXT,
        p_value REAL,
        odds_ratio REAL,
        beta REAL,
        ci_lower REAL,
        ci_upper REAL,
        risk_allele TEXT,
        gene TEXT,
        pmid TEXT,
        study_sample_size INTEGER,
        ancestry TEXT
    );
    
    CREATE INDEX IF NOT EXISTS idx_gwas_rsid ON associations(rsid);
    CREATE INDEX IF NOT EXISTS idx_gwas_trait ON associations(trait);
    CREATE INDEX IF NOT EXISTS idx_gwas_gene ON associations(gene);
    """
    
    def download(self, force: bool = False) -> bool:
        """Initialize GWAS Catalog with curated associations."""
        if self.is_downloaded and not force:
            logger.info("GWAS Catalog already initialized")
            return True
        
        logger.info("Initializing GWAS Catalog...")
        
        try:
            conn = self._get_connection()
            conn.executescript(self.SCHEMA)
            conn.commit()
            
            self._load_associations()
            
            version_info = DatasetVersion(
                name=self.name,
                version=self.version,
                downloaded=datetime.now(),
                source_url=self.source_url,
                record_count=self._count_associations(),
            )
            self.save_version_info(version_info)
            
            logger.info(f"GWAS Catalog initialized with {version_info.record_count} associations")
            return True
            
        except Exception as e:
            logger.error(f"Failed to initialize GWAS Catalog: {e}")
            return False
    
    def _load_associations(self) -> None:
        """Load curated GWAS associations."""
        ASSOCIATIONS = [
            # Height
            {"rsid": "rs1042725", "trait": "Height", "ontology": "EFO:0004339",
             "p_value": 1e-45, "beta": 0.35, "allele": "C", "gene": "HMGA2",
             "pmid": "25282103", "n": 253288, "ancestry": "European"},
            {"rsid": "rs3791679", "trait": "Height", "ontology": "EFO:0004339",
             "p_value": 3e-30, "beta": 0.28, "allele": "A", "gene": "LCORL",
             "pmid": "25282103", "n": 253288, "ancestry": "European"},
            
            # BMI/Obesity
            {"rsid": "rs1558902", "trait": "Body mass index", "ontology": "EFO:0004340",
             "p_value": 4e-120, "beta": 0.39, "allele": "A", "gene": "FTO",
             "pmid": "25673413", "n": 339224, "ancestry": "European"},
            {"rsid": "rs13107325", "trait": "Body mass index", "ontology": "EFO:0004340",
             "p_value": 2e-22, "beta": 0.19, "allele": "T", "gene": "SLC39A8",
             "pmid": "25673413", "n": 339224, "ancestry": "European"},
            
            # Type 2 Diabetes
            {"rsid": "rs7903146", "trait": "Type 2 diabetes", "ontology": "EFO:0001360",
             "p_value": 1e-350, "or": 1.40, "allele": "T", "gene": "TCF7L2",
             "pmid": "17463249", "n": 90000, "ancestry": "European"},
            {"rsid": "rs5219", "trait": "Type 2 diabetes", "ontology": "EFO:0001360",
             "p_value": 1e-10, "or": 1.15, "allele": "T", "gene": "KCNJ11",
             "pmid": "17463249", "n": 90000, "ancestry": "European"},
            
            # Coronary artery disease
            {"rsid": "rs1333049", "trait": "Coronary artery disease", "ontology": "EFO:0001645",
             "p_value": 1e-20, "or": 1.28, "allele": "C", "gene": "CDKN2A/B",
             "pmid": "17478679", "n": 23000, "ancestry": "European"},
            {"rsid": "rs10455872", "trait": "Coronary artery disease", "ontology": "EFO:0001645",
             "p_value": 1e-25, "or": 1.51, "allele": "G", "gene": "LPA",
             "pmid": "19060906", "n": 40000, "ancestry": "European"},
            
            # Alzheimer's disease
            {"rsid": "rs429358", "trait": "Alzheimer disease", "ontology": "EFO:0000249",
             "p_value": 1e-500, "or": 3.68, "allele": "C", "gene": "APOE",
             "pmid": "8346443", "n": 74000, "ancestry": "European"},
            {"rsid": "rs6656401", "trait": "Alzheimer disease", "ontology": "EFO:0000249",
             "p_value": 1e-16, "or": 1.18, "allele": "A", "gene": "CR1",
             "pmid": "21460841", "n": 56000, "ancestry": "European"},
            
            # Rheumatoid arthritis
            {"rsid": "rs2476601", "trait": "Rheumatoid arthritis", "ontology": "EFO:0000685",
             "p_value": 1e-35, "or": 1.75, "allele": "A", "gene": "PTPN22",
             "pmid": "15208781", "n": 20000, "ancestry": "European"},
            
            # Celiac disease
            {"rsid": "rs2187668", "trait": "Celiac disease", "ontology": "EFO:0001060",
             "p_value": 1e-200, "or": 7.0, "allele": "T", "gene": "HLA-DQA1",
             "pmid": "20190752", "n": 15283, "ancestry": "European"},
            
            # Multiple sclerosis
            {"rsid": "rs3135388", "trait": "Multiple sclerosis", "ontology": "EFO:0003885",
             "p_value": 1e-100, "or": 2.5, "allele": "A", "gene": "HLA-DRB1",
             "pmid": "17660530", "n": 30000, "ancestry": "European"},
            
            # Schizophrenia
            {"rsid": "rs2021722", "trait": "Schizophrenia", "ontology": "EFO:0000692",
             "p_value": 1e-20, "or": 1.12, "allele": "T", "gene": "C4A",
             "pmid": "25056061", "n": 150000, "ancestry": "European"},
            
            # Bipolar disorder
            {"rsid": "rs4765913", "trait": "Bipolar disorder", "ontology": "EFO:0000289",
             "p_value": 1e-12, "or": 1.12, "allele": "C", "gene": "CACNA1C",
             "pmid": "30804558", "n": 30000, "ancestry": "European"},
            
            # Age-related macular degeneration
            {"rsid": "rs1061170", "trait": "Age-related macular degeneration", "ontology": "EFO:0001365",
             "p_value": 1e-100, "or": 2.5, "allele": "C", "gene": "CFH",
             "pmid": "15761122", "n": 30000, "ancestry": "European"},
            
            # Caffeine metabolism
            {"rsid": "rs762551", "trait": "Caffeine metabolism", "ontology": "EFO:0007797",
             "p_value": 1e-50, "beta": -0.5, "allele": "C", "gene": "CYP1A2",
             "pmid": "22426360", "n": 120000, "ancestry": "European"},
            
            # Coffee consumption
            {"rsid": "rs2472297", "trait": "Coffee consumption", "ontology": "EFO:0007796",
             "p_value": 1e-60, "beta": 0.15, "allele": "T", "gene": "CYP1A2",
             "pmid": "25288136", "n": 120000, "ancestry": "European"},
            
            # Alcohol consumption
            {"rsid": "rs1229984", "trait": "Alcohol consumption", "ontology": "EFO:0007878",
             "p_value": 1e-100, "beta": -0.5, "allele": "A", "gene": "ADH1B",
             "pmid": "31358974", "n": 480000, "ancestry": "Multiple"},
            
            # Sleep duration
            {"rsid": "rs1805247", "trait": "Sleep duration", "ontology": "EFO:0005271",
             "p_value": 1e-20, "beta": 0.05, "allele": "T", "gene": "GABRA2",
             "pmid": "27618452", "n": 450000, "ancestry": "European"},
            
            # Chronotype (morning/evening preference)
            {"rsid": "rs12736689", "trait": "Chronotype", "ontology": "EFO:0004574",
             "p_value": 1e-40, "beta": 0.10, "allele": "G", "gene": "PER2",
             "pmid": "27494321", "n": 450000, "ancestry": "European"},
            
            # Eye color
            {"rsid": "rs12913832", "trait": "Eye color", "ontology": "EFO:0003927",
             "p_value": 1e-300, "or": 30.0, "allele": "G", "gene": "HERC2",
             "pmid": "18488028", "n": 6000, "ancestry": "European"},
            
            # Hair color
            {"rsid": "rs1805007", "trait": "Hair color", "ontology": "EFO:0003926",
             "p_value": 1e-50, "or": 5.0, "allele": "T", "gene": "MC1R",
             "pmid": "18488028", "n": 6000, "ancestry": "European"},
            
            # Lactose intolerance
            {"rsid": "rs4988235", "trait": "Lactose intolerance", "ontology": "EFO:0004272",
             "p_value": 1e-200, "or": 0.05, "allele": "G", "gene": "MCM6",
             "pmid": "12858176", "n": 10000, "ancestry": "European"},
            
            # Asparagus anosmia (ability to smell)
            {"rsid": "rs4481887", "trait": "Asparagus anosmia", "ontology": "EFO:0008330",
             "p_value": 1e-30, "or": 2.5, "allele": "G", "gene": "OR2M7",
             "pmid": "27911795", "n": 9000, "ancestry": "European"},
            
            # Bitter taste perception
            {"rsid": "rs713598", "trait": "Bitter taste perception", "ontology": "EFO:0004252",
             "p_value": 1e-50, "or": 5.0, "allele": "C", "gene": "TAS2R38",
             "pmid": "12595690", "n": 5000, "ancestry": "European"},
            
            # Cilantro soapy taste
            {"rsid": "rs72921001", "trait": "Cilantro preference", "ontology": "EFO:0008334",
             "p_value": 1e-15, "or": 0.4, "allele": "A", "gene": "OR6A2",
             "pmid": "22927436", "n": 30000, "ancestry": "European"},
        ]
        
        conn = self._get_connection()
        cursor = conn.cursor()
        
        for assoc in ASSOCIATIONS:
            cursor.execute("""
                INSERT INTO associations 
                (rsid, trait, trait_ontology, p_value, odds_ratio, beta, 
                 risk_allele, gene, pmid, study_sample_size, ancestry)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (assoc["rsid"], assoc["trait"], assoc["ontology"],
                  assoc["p_value"], assoc.get("or"), assoc.get("beta"),
                  assoc["allele"], assoc["gene"], assoc["pmid"],
                  assoc["n"], assoc["ancestry"]))
        
        conn.commit()
        logger.info(f"Loaded {len(ASSOCIATIONS)} GWAS associations")
    
    def _count_associations(self) -> int:
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM associations")
        return cursor.fetchone()[0]
    
    def lookup_variant(self, rsid: str) -> Optional[VariantInfo]:
        """Look up GWAS associations for a variant."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT DISTINCT gene FROM associations WHERE rsid = ?
        """, (rsid,))
        row = cursor.fetchone()
        
        if not row:
            return None
        
        return VariantInfo(
            rsid=rsid,
            chromosome="",
            position=0,
            ref_allele="",
            alt_allele="",
            gene=row["gene"],
        )
    
    def get_variant_associations(self, rsid: str) -> List[GWASAssociation]:
        """Get all GWAS associations for a variant."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("SELECT * FROM associations WHERE rsid = ?", (rsid,))
        
        results = []
        for row in cursor.fetchall():
            results.append(GWASAssociation(
                rsid=row["rsid"],
                trait=row["trait"],
                trait_ontology=row["trait_ontology"],
                p_value=row["p_value"],
                odds_ratio=row["odds_ratio"],
                beta=row["beta"],
                ci_lower=row["ci_lower"],
                ci_upper=row["ci_upper"],
                risk_allele=row["risk_allele"],
                gene=row["gene"],
                pmid=row["pmid"],
                study_sample_size=row["study_sample_size"],
                ancestry=row["ancestry"],
            ))
        
        return results
    
    def get_trait_variants(self, trait: str) -> List[GWASAssociation]:
        """Get all variants associated with a trait."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT * FROM associations WHERE LOWER(trait) LIKE ?
        """, (f"%{trait.lower()}%",))
        
        results = []
        for row in cursor.fetchall():
            results.append(GWASAssociation(
                rsid=row["rsid"],
                trait=row["trait"],
                trait_ontology=row["trait_ontology"],
                p_value=row["p_value"],
                odds_ratio=row["odds_ratio"],
                beta=row["beta"],
                ci_lower=row["ci_lower"],
                ci_upper=row["ci_upper"],
                risk_allele=row["risk_allele"],
                gene=row["gene"],
                pmid=row["pmid"],
                study_sample_size=row["study_sample_size"],
                ancestry=row["ancestry"],
            ))
        
        return results
    
    def check_user_variants(
        self, 
        genotypes: Dict[str, str]
    ) -> Dict[str, List[GWASAssociation]]:
        """
        Check user's variants against GWAS Catalog.
        
        Returns dict mapping rsid to list of associations.
        """
        results = {}
        
        for rsid in genotypes:
            associations = self.get_variant_associations(rsid)
            if associations:
                results[rsid] = associations
        
        return results


# =============================================================================
# CONVENIENCE FUNCTIONS  
# =============================================================================

def get_trait_associations(rsid: str) -> List[GWASAssociation]:
    """Get GWAS associations for a variant."""
    gwas = GWASCatalog()
    if not gwas.is_downloaded:
        gwas.download()
    return gwas.get_variant_associations(rsid)
