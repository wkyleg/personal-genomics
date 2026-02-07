"""
PharmGKB (Pharmacogenomics Knowledge Base) Integration

Provides drug-gene interaction annotations and dosing guidelines.
Critical for medication safety and pharmacogenomics analysis.

Data source: https://www.pharmgkb.org/
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


# PharmGKB evidence levels
EVIDENCE_LEVELS = {
    "1A": "Strong, implemented in a clinical guideline",
    "1B": "Strong, clinical guideline annotation",
    "2A": "Moderate, functional significance",
    "2B": "Moderate, replicated associations",
    "3": "Low, single study or conflicting evidence",
    "4": "Preliminary, case reports or in vitro only",
}

# CPIC (Clinical Pharmacogenetics Implementation Consortium) guidelines
# These are the gold standard for pharmacogenomics
CPIC_GUIDELINES = {
    "CYP2D6": ["codeine", "tramadol", "oxycodone", "hydrocodone", "tamoxifen", 
               "ondansetron", "tropisetron", "atomoxetine", "venlafaxine",
               "desipramine", "nortriptyline", "amitriptyline", "clomipramine"],
    "CYP2C19": ["clopidogrel", "voriconazole", "proton pump inhibitors", 
                "amitriptyline", "citalopram", "escitalopram", "sertraline"],
    "CYP2C9": ["warfarin", "phenytoin"],
    "VKORC1": ["warfarin"],
    "TPMT": ["azathioprine", "mercaptopurine", "thioguanine"],
    "NUDT15": ["azathioprine", "mercaptopurine", "thioguanine"],
    "DPYD": ["fluorouracil", "capecitabine", "tegafur"],
    "UGT1A1": ["irinotecan", "atazanavir"],
    "SLCO1B1": ["simvastatin", "atorvastatin", "rosuvastatin", "pravastatin"],
    "CYP3A5": ["tacrolimus"],
    "HLA-B": ["abacavir", "carbamazepine", "oxcarbazepine", "phenytoin", "allopurinol"],
    "HLA-A": ["carbamazepine"],
    "G6PD": ["rasburicase", "dapsone", "primaquine"],
    "CYP4F2": ["warfarin"],
    "IFNL3": ["peginterferon alfa-2a", "peginterferon alfa-2b"],
    "RYR1": ["volatile anesthetics", "succinylcholine"],
    "CACNA1S": ["volatile anesthetics", "succinylcholine"],
}


@dataclass
class DrugGeneInteraction:
    """PharmGKB drug-gene interaction annotation."""
    gene: str
    drug: str
    rsid: Optional[str]
    variant_name: Optional[str]  # e.g., CYP2D6*4
    phenotype: str  # e.g., "Poor Metabolizer"
    recommendation: str
    evidence_level: str
    cpic_level: Optional[str]
    pmids: List[str] = field(default_factory=list)
    
    @property
    def evidence_description(self) -> str:
        return EVIDENCE_LEVELS.get(self.evidence_level, "Unknown")
    
    @property
    def is_actionable(self) -> bool:
        """Check if this interaction has clinical actionability."""
        return self.evidence_level in ["1A", "1B", "2A"]


@dataclass
class DosingGuideline:
    """CPIC dosing guideline."""
    gene: str
    drug: str
    phenotype: str
    recommendation: str
    alternative_drugs: List[str]
    dose_adjustment: Optional[str]
    cpic_level: str
    strength: str  # "strong", "moderate", "optional"


class PharmGKB(SQLiteDataset):
    """
    PharmGKB pharmacogenomics database integration.
    
    Provides:
    - Drug-gene interactions with evidence levels
    - CPIC dosing guidelines
    - Star allele definitions
    - Phenotype predictions
    """
    
    name = "pharmgkb"
    version = "2024-02"
    description = "PharmGKB pharmacogenomics annotations"
    source_url = "https://www.pharmgkb.org/"
    
    SCHEMA = """
    CREATE TABLE IF NOT EXISTS drug_gene_interactions (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene TEXT,
        drug TEXT,
        rsid TEXT,
        variant_name TEXT,
        phenotype TEXT,
        recommendation TEXT,
        evidence_level TEXT,
        cpic_level TEXT
    );
    
    CREATE TABLE IF NOT EXISTS dosing_guidelines (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene TEXT,
        drug TEXT,
        phenotype TEXT,
        recommendation TEXT,
        dose_adjustment TEXT,
        cpic_level TEXT,
        strength TEXT
    );
    
    CREATE TABLE IF NOT EXISTS alternative_drugs (
        guideline_id INTEGER,
        alternative_drug TEXT,
        FOREIGN KEY (guideline_id) REFERENCES dosing_guidelines(id)
    );
    
    CREATE TABLE IF NOT EXISTS variant_definitions (
        rsid TEXT PRIMARY KEY,
        gene TEXT,
        star_allele TEXT,
        function TEXT,
        activity_value REAL
    );
    
    CREATE TABLE IF NOT EXISTS interaction_pmids (
        interaction_id INTEGER,
        pmid TEXT,
        FOREIGN KEY (interaction_id) REFERENCES drug_gene_interactions(id)
    );
    
    CREATE INDEX IF NOT EXISTS idx_pharmgkb_gene ON drug_gene_interactions(gene);
    CREATE INDEX IF NOT EXISTS idx_pharmgkb_drug ON drug_gene_interactions(drug);
    CREATE INDEX IF NOT EXISTS idx_pharmgkb_rsid ON drug_gene_interactions(rsid);
    CREATE INDEX IF NOT EXISTS idx_dosing_gene ON dosing_guidelines(gene);
    CREATE INDEX IF NOT EXISTS idx_dosing_drug ON dosing_guidelines(drug);
    """
    
    def download(self, force: bool = False) -> bool:
        """Initialize PharmGKB database with CPIC guidelines."""
        if self.is_downloaded and not force:
            logger.info("PharmGKB data already downloaded")
            return True
        
        logger.info("Initializing PharmGKB data...")
        
        try:
            conn = self._get_connection()
            conn.executescript(self.SCHEMA)
            conn.commit()
            
            # Load CPIC guidelines and variant definitions
            self._load_cpic_guidelines()
            self._load_variant_definitions()
            
            version_info = DatasetVersion(
                name=self.name,
                version=self.version,
                downloaded=datetime.now(),
                source_url=self.source_url,
                record_count=self._count_interactions(),
            )
            self.save_version_info(version_info)
            
            logger.info(f"PharmGKB initialized with {version_info.record_count} interactions")
            return True
            
        except Exception as e:
            logger.error(f"Failed to initialize PharmGKB: {e}")
            return False
    
    def _load_cpic_guidelines(self) -> None:
        """Load CPIC dosing guidelines."""
        GUIDELINES = [
            # CYP2D6 - Codeine
            {
                "gene": "CYP2D6", "drug": "codeine",
                "phenotypes": {
                    "Ultrarapid Metabolizer": {
                        "recommendation": "AVOID codeine. Use alternative analgesic (not tramadol/oxycodone)",
                        "dose_adjustment": "Do not use",
                        "strength": "strong",
                        "alternatives": ["morphine", "non-opioid analgesics", "NSAIDs"],
                    },
                    "Poor Metabolizer": {
                        "recommendation": "AVOID codeine due to lack of efficacy. Use alternative analgesic",
                        "dose_adjustment": "Do not use",
                        "strength": "strong",
                        "alternatives": ["morphine", "non-opioid analgesics", "NSAIDs"],
                    },
                    "Intermediate Metabolizer": {
                        "recommendation": "Use codeine with caution or consider alternative",
                        "dose_adjustment": "Reduce dose by 25-50%",
                        "strength": "moderate",
                        "alternatives": ["morphine", "non-opioid analgesics"],
                    },
                    "Normal Metabolizer": {
                        "recommendation": "Use codeine per standard dosing",
                        "dose_adjustment": None,
                        "strength": "optional",
                        "alternatives": [],
                    },
                }
            },
            # CYP2C19 - Clopidogrel
            {
                "gene": "CYP2C19", "drug": "clopidogrel",
                "phenotypes": {
                    "Poor Metabolizer": {
                        "recommendation": "Use alternative antiplatelet (prasugrel/ticagrelor if no contraindications)",
                        "dose_adjustment": "Do not use",
                        "strength": "strong",
                        "alternatives": ["prasugrel", "ticagrelor"],
                    },
                    "Intermediate Metabolizer": {
                        "recommendation": "Consider alternative antiplatelet, especially for high-risk conditions (ACS, PCI)",
                        "dose_adjustment": "Consider alternative",
                        "strength": "moderate",
                        "alternatives": ["prasugrel", "ticagrelor"],
                    },
                    "Normal Metabolizer": {
                        "recommendation": "Use clopidogrel per standard dosing",
                        "dose_adjustment": None,
                        "strength": "optional",
                        "alternatives": [],
                    },
                    "Rapid Metabolizer": {
                        "recommendation": "Use clopidogrel per standard dosing",
                        "dose_adjustment": None,
                        "strength": "optional",
                        "alternatives": [],
                    },
                }
            },
            # DPYD - Fluorouracil
            {
                "gene": "DPYD", "drug": "fluorouracil",
                "phenotypes": {
                    "Poor Metabolizer": {
                        "recommendation": "AVOID fluorouracil/capecitabine - HIGH TOXICITY RISK (potentially fatal)",
                        "dose_adjustment": "Do not use",
                        "strength": "strong",
                        "alternatives": ["alternative chemotherapy regimen"],
                    },
                    "Intermediate Metabolizer": {
                        "recommendation": "Reduce dose by 50% with intensive monitoring",
                        "dose_adjustment": "Reduce by 50%",
                        "strength": "strong",
                        "alternatives": [],
                    },
                    "Normal Metabolizer": {
                        "recommendation": "Use standard dosing",
                        "dose_adjustment": None,
                        "strength": "optional",
                        "alternatives": [],
                    },
                }
            },
            # SLCO1B1 - Simvastatin
            {
                "gene": "SLCO1B1", "drug": "simvastatin",
                "phenotypes": {
                    "Poor Function": {
                        "recommendation": "Use lower dose simvastatin (20mg max) or alternative statin",
                        "dose_adjustment": "Maximum 20mg daily",
                        "strength": "strong",
                        "alternatives": ["pravastatin", "rosuvastatin"],
                    },
                    "Intermediate Function": {
                        "recommendation": "Consider lower dose or alternative if myopathy risk factors present",
                        "dose_adjustment": "Consider reduction",
                        "strength": "moderate",
                        "alternatives": ["pravastatin", "rosuvastatin"],
                    },
                    "Normal Function": {
                        "recommendation": "Use standard dosing",
                        "dose_adjustment": None,
                        "strength": "optional",
                        "alternatives": [],
                    },
                }
            },
            # HLA-B*5701 - Abacavir
            {
                "gene": "HLA-B", "drug": "abacavir",
                "phenotypes": {
                    "HLA-B*5701 Positive": {
                        "recommendation": "DO NOT USE abacavir - hypersensitivity reaction risk",
                        "dose_adjustment": "Contraindicated",
                        "strength": "strong",
                        "alternatives": ["other NRTIs", "tenofovir"],
                    },
                    "HLA-B*5701 Negative": {
                        "recommendation": "Use abacavir per standard dosing",
                        "dose_adjustment": None,
                        "strength": "optional",
                        "alternatives": [],
                    },
                }
            },
            # TPMT - Azathioprine/Mercaptopurine
            {
                "gene": "TPMT", "drug": "azathioprine",
                "phenotypes": {
                    "Poor Metabolizer": {
                        "recommendation": "Start with drastically reduced dose (10% of normal) or avoid",
                        "dose_adjustment": "10% of normal dose",
                        "strength": "strong",
                        "alternatives": [],
                    },
                    "Intermediate Metabolizer": {
                        "recommendation": "Start with reduced dose (30-70% of normal)",
                        "dose_adjustment": "30-70% of normal",
                        "strength": "strong",
                        "alternatives": [],
                    },
                    "Normal Metabolizer": {
                        "recommendation": "Use standard dosing",
                        "dose_adjustment": None,
                        "strength": "optional",
                        "alternatives": [],
                    },
                }
            },
            # CYP2C9/VKORC1 - Warfarin
            {
                "gene": "CYP2C9", "drug": "warfarin",
                "phenotypes": {
                    "Poor Metabolizer": {
                        "recommendation": "Significantly reduced dose required (use FDA table or IWPC algorithm)",
                        "dose_adjustment": "Reduce by 50-80%",
                        "strength": "strong",
                        "alternatives": ["direct oral anticoagulants"],
                    },
                    "Intermediate Metabolizer": {
                        "recommendation": "Reduced dose likely required",
                        "dose_adjustment": "Reduce by 20-40%",
                        "strength": "moderate",
                        "alternatives": [],
                    },
                    "Normal Metabolizer": {
                        "recommendation": "Use standard dosing algorithm",
                        "dose_adjustment": None,
                        "strength": "optional",
                        "alternatives": [],
                    },
                }
            },
        ]
        
        conn = self._get_connection()
        cursor = conn.cursor()
        
        for guideline in GUIDELINES:
            gene = guideline["gene"]
            drug = guideline["drug"]
            
            for phenotype, details in guideline["phenotypes"].items():
                cursor.execute("""
                    INSERT INTO dosing_guidelines 
                    (gene, drug, phenotype, recommendation, dose_adjustment, cpic_level, strength)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """, (gene, drug, phenotype, details["recommendation"],
                      details["dose_adjustment"], "1A", details["strength"]))
                
                guideline_id = cursor.lastrowid
                
                for alt_drug in details.get("alternatives", []):
                    cursor.execute("""
                        INSERT INTO alternative_drugs (guideline_id, alternative_drug)
                        VALUES (?, ?)
                    """, (guideline_id, alt_drug))
        
        conn.commit()
        logger.info(f"Loaded {len(GUIDELINES)} CPIC guidelines")
    
    def _load_variant_definitions(self) -> None:
        """Load pharmacogene variant definitions (star alleles)."""
        VARIANTS = {
            # CYP2D6 star alleles
            "rs3892097": {"gene": "CYP2D6", "star": "*4", "function": "No function", "activity": 0.0},
            "rs5030655": {"gene": "CYP2D6", "star": "*6", "function": "No function", "activity": 0.0},
            "rs16947": {"gene": "CYP2D6", "star": "*2", "function": "Normal function", "activity": 1.0},
            "rs1065852": {"gene": "CYP2D6", "star": "*10", "function": "Decreased function", "activity": 0.25},
            "rs28371706": {"gene": "CYP2D6", "star": "*17", "function": "Decreased function", "activity": 0.5},
            "rs1135840": {"gene": "CYP2D6", "star": "*2/*41", "function": "Normal/Decreased", "activity": 0.75},
            
            # CYP2C19 star alleles
            "rs4244285": {"gene": "CYP2C19", "star": "*2", "function": "No function", "activity": 0.0},
            "rs4986893": {"gene": "CYP2C19", "star": "*3", "function": "No function", "activity": 0.0},
            "rs12248560": {"gene": "CYP2C19", "star": "*17", "function": "Increased function", "activity": 1.5},
            
            # CYP2C9 star alleles
            "rs1799853": {"gene": "CYP2C9", "star": "*2", "function": "Decreased function", "activity": 0.5},
            "rs1057910": {"gene": "CYP2C9", "star": "*3", "function": "Decreased function", "activity": 0.1},
            
            # DPYD variants
            "rs3918290": {"gene": "DPYD", "star": "*2A", "function": "No function", "activity": 0.0},
            "rs55886062": {"gene": "DPYD", "star": "c.1679T>G", "function": "Decreased function", "activity": 0.0},
            "rs67376798": {"gene": "DPYD", "star": "D949V", "function": "Decreased function", "activity": 0.5},
            
            # TPMT variants
            "rs1800460": {"gene": "TPMT", "star": "*3A/*3C", "function": "No function", "activity": 0.0},
            "rs1142345": {"gene": "TPMT", "star": "*3A/*3B", "function": "No function", "activity": 0.0},
            
            # SLCO1B1 variants
            "rs4149056": {"gene": "SLCO1B1", "star": "*5", "function": "Decreased function", "activity": 0.3},
            "rs2306283": {"gene": "SLCO1B1", "star": "*1B", "function": "Normal function", "activity": 1.0},
            
            # VKORC1 variants
            "rs9923231": {"gene": "VKORC1", "star": "-1639G>A", "function": "Decreased expression", "activity": 0.5},
            
            # UGT1A1 variants
            "rs8175347": {"gene": "UGT1A1", "star": "*28", "function": "Decreased function", "activity": 0.3},
        }
        
        conn = self._get_connection()
        cursor = conn.cursor()
        
        for rsid, data in VARIANTS.items():
            cursor.execute("""
                INSERT OR REPLACE INTO variant_definitions 
                (rsid, gene, star_allele, function, activity_value)
                VALUES (?, ?, ?, ?, ?)
            """, (rsid, data["gene"], data["star"], data["function"], data["activity"]))
        
        conn.commit()
        logger.info(f"Loaded {len(VARIANTS)} pharmacogene variants")
    
    def _count_interactions(self) -> int:
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM dosing_guidelines")
        return cursor.fetchone()[0]
    
    def lookup_variant(self, rsid: str) -> Optional[VariantInfo]:
        """Look up variant definition by rsID."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT rsid, gene, star_allele, function, activity_value
            FROM variant_definitions WHERE rsid = ?
        """, (rsid,))
        row = cursor.fetchone()
        
        if not row:
            return None
        
        return VariantInfo(
            rsid=row["rsid"],
            chromosome="",
            position=0,
            ref_allele="",
            alt_allele="",
            gene=row["gene"],
        )
    
    def get_drug_interactions(self, gene: str) -> List[DosingGuideline]:
        """Get all drug interactions for a gene."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT id, gene, drug, phenotype, recommendation, dose_adjustment, 
                   cpic_level, strength
            FROM dosing_guidelines WHERE gene = ?
        """, (gene,))
        
        results = []
        for row in cursor.fetchall():
            cursor.execute("""
                SELECT alternative_drug FROM alternative_drugs WHERE guideline_id = ?
            """, (row["id"],))
            alternatives = [r["alternative_drug"] for r in cursor.fetchall()]
            
            results.append(DosingGuideline(
                gene=row["gene"],
                drug=row["drug"],
                phenotype=row["phenotype"],
                recommendation=row["recommendation"],
                alternative_drugs=alternatives,
                dose_adjustment=row["dose_adjustment"],
                cpic_level=row["cpic_level"],
                strength=row["strength"],
            ))
        
        return results
    
    def get_dosing_recommendation(
        self, 
        gene: str, 
        drug: str, 
        phenotype: str
    ) -> Optional[DosingGuideline]:
        """Get specific dosing recommendation."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT id, gene, drug, phenotype, recommendation, dose_adjustment, 
                   cpic_level, strength
            FROM dosing_guidelines 
            WHERE gene = ? AND drug = ? AND phenotype = ?
        """, (gene, drug, phenotype))
        
        row = cursor.fetchone()
        if not row:
            return None
        
        cursor.execute("""
            SELECT alternative_drug FROM alternative_drugs WHERE guideline_id = ?
        """, (row["id"],))
        alternatives = [r["alternative_drug"] for r in cursor.fetchall()]
        
        return DosingGuideline(
            gene=row["gene"],
            drug=row["drug"],
            phenotype=row["phenotype"],
            recommendation=row["recommendation"],
            alternative_drugs=alternatives,
            dose_adjustment=row["dose_adjustment"],
            cpic_level=row["cpic_level"],
            strength=row["strength"],
        )
    
    def get_star_allele(self, rsid: str) -> Optional[Tuple[str, str, float]]:
        """Get star allele assignment for a variant."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT gene, star_allele, activity_value
            FROM variant_definitions WHERE rsid = ?
        """, (rsid,))
        row = cursor.fetchone()
        
        if not row:
            return None
        
        return (row["gene"], row["star_allele"], row["activity_value"])
    
    def calculate_activity_score(
        self, 
        gene: str, 
        genotypes: Dict[str, str]
    ) -> Tuple[float, str]:
        """
        Calculate activity score for a pharmacogene based on genotypes.
        
        Returns (activity_score, phenotype) tuple.
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        
        # Get all variants for this gene
        cursor.execute("""
            SELECT rsid, star_allele, function, activity_value
            FROM variant_definitions WHERE gene = ?
        """, (gene,))
        
        gene_variants = {row["rsid"]: row for row in cursor.fetchall()}
        
        if not gene_variants:
            return (2.0, "Normal Metabolizer")  # Default assumption
        
        total_activity = 0.0
        variants_found = 0
        
        for rsid, var_info in gene_variants.items():
            geno = genotypes.get(rsid)
            if not geno:
                continue
            
            variants_found += 1
            activity = var_info["activity_value"]
            
            # Count variant alleles (simplified diplotype assignment)
            # Real implementation would need proper haplotype phasing
            # For now, assume heterozygous means one copy of variant
            if geno[0] != geno[1]:  # Heterozygous
                total_activity += 1.0 + activity  # One normal + one variant
            else:
                total_activity += 2 * activity  # Homozygous variant
        
        if variants_found == 0:
            return (2.0, "Normal Metabolizer")  # Assume normal if no data
        
        # Normalize to diploid activity score
        avg_activity = total_activity / variants_found
        
        # Assign phenotype based on activity score
        if avg_activity >= 2.5:
            phenotype = "Ultrarapid Metabolizer"
        elif avg_activity >= 1.5:
            phenotype = "Rapid Metabolizer"
        elif avg_activity >= 1.0:
            phenotype = "Normal Metabolizer"
        elif avg_activity >= 0.5:
            phenotype = "Intermediate Metabolizer"
        else:
            phenotype = "Poor Metabolizer"
        
        return (avg_activity, phenotype)


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def get_drug_recommendations(gene: str) -> List[DosingGuideline]:
    """Get CPIC dosing recommendations for a gene."""
    pgkb = PharmGKB()
    if not pgkb.is_downloaded:
        pgkb.download()
    return pgkb.get_drug_interactions(gene)


def check_medication_safety(
    medications: List[str], 
    genotypes: Dict[str, str]
) -> Dict[str, List[DosingGuideline]]:
    """Check medications against pharmacogenomics guidelines."""
    pgkb = PharmGKB()
    if not pgkb.is_downloaded:
        pgkb.download()
    
    results = {}
    
    for gene in CPIC_GUIDELINES:
        activity_score, phenotype = pgkb.calculate_activity_score(gene, genotypes)
        
        for med in medications:
            med_lower = med.lower()
            if any(drug.lower() in med_lower or med_lower in drug.lower() 
                   for drug in CPIC_GUIDELINES.get(gene, [])):
                rec = pgkb.get_dosing_recommendation(gene, med_lower, phenotype)
                if rec:
                    if med not in results:
                        results[med] = []
                    results[med].append(rec)
    
    return results
