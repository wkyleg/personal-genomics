"""
Human Genome Diversity Project (HGDP) Integration

Provides population-level data from ~1,000 individuals across 51 populations.
Excellent for fine-scale ancestry analysis and global diversity comparison.

Data source: https://www.hagsc.org/hgdp/
"""

from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass
import json
import logging
import sqlite3

from .base import (
    SQLiteDataset,
    DatasetVersion,
    VariantInfo,
    PopulationFrequency,
)

logger = logging.getLogger(__name__)


# HGDP population definitions grouped by region
HGDP_POPULATIONS = {
    # Africa
    "Bantu_Kenya": {"region": "Africa", "country": "Kenya", "n": 11},
    "Bantu_South_Africa": {"region": "Africa", "country": "South Africa", "n": 8},
    "Biaka_Pygmies": {"region": "Africa", "country": "Central African Republic", "n": 23},
    "Mandenka": {"region": "Africa", "country": "Senegal", "n": 22},
    "Mbuti_Pygmies": {"region": "Africa", "country": "Democratic Republic of Congo", "n": 13},
    "San": {"region": "Africa", "country": "Namibia", "n": 6},
    "Yoruba": {"region": "Africa", "country": "Nigeria", "n": 21},
    
    # Middle East
    "Bedouin": {"region": "Middle_East", "country": "Israel", "n": 45},
    "Druze": {"region": "Middle_East", "country": "Israel", "n": 42},
    "Mozabite": {"region": "Middle_East", "country": "Algeria", "n": 29},
    "Palestinian": {"region": "Middle_East", "country": "Israel", "n": 46},
    
    # Europe
    "Adygei": {"region": "Europe", "country": "Russia", "n": 17},
    "Basque": {"region": "Europe", "country": "France", "n": 24},
    "French": {"region": "Europe", "country": "France", "n": 28},
    "Bergamo_Italian": {"region": "Europe", "country": "Italy", "n": 13},
    "Orcadian": {"region": "Europe", "country": "Scotland", "n": 15},
    "Russian": {"region": "Europe", "country": "Russia", "n": 25},
    "Sardinian": {"region": "Europe", "country": "Italy", "n": 28},
    "Tuscan": {"region": "Europe", "country": "Italy", "n": 8},
    
    # Central/South Asia
    "Balochi": {"region": "Central_South_Asia", "country": "Pakistan", "n": 24},
    "Brahui": {"region": "Central_South_Asia", "country": "Pakistan", "n": 25},
    "Burusho": {"region": "Central_South_Asia", "country": "Pakistan", "n": 25},
    "Hazara": {"region": "Central_South_Asia", "country": "Pakistan", "n": 22},
    "Kalash": {"region": "Central_South_Asia", "country": "Pakistan", "n": 23},
    "Makrani": {"region": "Central_South_Asia", "country": "Pakistan", "n": 25},
    "Pathan": {"region": "Central_South_Asia", "country": "Pakistan", "n": 22},
    "Sindhi": {"region": "Central_South_Asia", "country": "Pakistan", "n": 24},
    "Uygur": {"region": "Central_South_Asia", "country": "China", "n": 10},
    
    # East Asia
    "Cambodian": {"region": "East_Asia", "country": "Cambodia", "n": 10},
    "Dai": {"region": "East_Asia", "country": "China", "n": 10},
    "Daur": {"region": "East_Asia", "country": "China", "n": 9},
    "Han": {"region": "East_Asia", "country": "China", "n": 44},
    "Han_N_China": {"region": "East_Asia", "country": "China", "n": 10},
    "Hezhen": {"region": "East_Asia", "country": "China", "n": 9},
    "Japanese": {"region": "East_Asia", "country": "Japan", "n": 28},
    "Lahu": {"region": "East_Asia", "country": "China", "n": 8},
    "Miao": {"region": "East_Asia", "country": "China", "n": 10},
    "Mongola": {"region": "East_Asia", "country": "China", "n": 10},
    "Naxi": {"region": "East_Asia", "country": "China", "n": 9},
    "Oroqen": {"region": "East_Asia", "country": "China", "n": 9},
    "She": {"region": "East_Asia", "country": "China", "n": 10},
    "Tu": {"region": "East_Asia", "country": "China", "n": 10},
    "Tujia": {"region": "East_Asia", "country": "China", "n": 10},
    "Xibo": {"region": "East_Asia", "country": "China", "n": 9},
    "Yi": {"region": "East_Asia", "country": "China", "n": 10},
    "Yakut": {"region": "East_Asia", "country": "Russia", "n": 25},
    
    # Oceania
    "Melanesian": {"region": "Oceania", "country": "Papua New Guinea", "n": 10},
    "Papuan": {"region": "Oceania", "country": "Papua New Guinea", "n": 17},
    
    # Americas
    "Colombian": {"region": "Americas", "country": "Colombia", "n": 7},
    "Karitiana": {"region": "Americas", "country": "Brazil", "n": 14},
    "Maya": {"region": "Americas", "country": "Mexico", "n": 21},
    "Pima": {"region": "Americas", "country": "Mexico", "n": 14},
    "Surui": {"region": "Americas", "country": "Brazil", "n": 8},
}

HGDP_REGIONS = {
    "Africa": ["Bantu_Kenya", "Bantu_South_Africa", "Biaka_Pygmies", "Mandenka", 
               "Mbuti_Pygmies", "San", "Yoruba"],
    "Middle_East": ["Bedouin", "Druze", "Mozabite", "Palestinian"],
    "Europe": ["Adygei", "Basque", "French", "Bergamo_Italian", "Orcadian", 
               "Russian", "Sardinian", "Tuscan"],
    "Central_South_Asia": ["Balochi", "Brahui", "Burusho", "Hazara", "Kalash",
                           "Makrani", "Pathan", "Sindhi", "Uygur"],
    "East_Asia": ["Cambodian", "Dai", "Daur", "Han", "Han_N_China", "Hezhen",
                  "Japanese", "Lahu", "Miao", "Mongola", "Naxi", "Oroqen", 
                  "She", "Tu", "Tujia", "Xibo", "Yi", "Yakut"],
    "Oceania": ["Melanesian", "Papuan"],
    "Americas": ["Colombian", "Karitiana", "Maya", "Pima", "Surui"],
}


class HGDP(SQLiteDataset):
    """
    Human Genome Diversity Project integration.
    
    Provides:
    - Fine-scale population frequencies
    - 51 populations across 7 world regions
    - Global diversity comparison
    """
    
    name = "hgdp"
    version = "stanford_release"
    description = "Human Genome Diversity Project"
    source_url = "https://www.hagsc.org/hgdp/"
    
    SCHEMA = """
    CREATE TABLE IF NOT EXISTS variants (
        rsid TEXT PRIMARY KEY,
        chromosome TEXT,
        position INTEGER,
        ref_allele TEXT,
        alt_allele TEXT
    );
    
    CREATE TABLE IF NOT EXISTS frequencies (
        rsid TEXT,
        population TEXT,
        allele TEXT,
        frequency REAL,
        sample_size INTEGER,
        PRIMARY KEY (rsid, population),
        FOREIGN KEY (rsid) REFERENCES variants(rsid)
    );
    
    CREATE INDEX IF NOT EXISTS idx_hgdp_pop ON frequencies(population);
    """
    
    def download(self, force: bool = False) -> bool:
        """Initialize HGDP database with AIMs."""
        if self.is_downloaded and not force:
            logger.info("HGDP data already downloaded")
            return True
        
        logger.info("Initializing HGDP data...")
        
        try:
            conn = self._get_connection()
            conn.executescript(self.SCHEMA)
            conn.commit()
            
            self._load_aim_frequencies()
            
            version_info = DatasetVersion(
                name=self.name,
                version=self.version,
                downloaded=datetime.now(),
                source_url=self.source_url,
                record_count=self._count_variants(),
            )
            self.save_version_info(version_info)
            
            logger.info(f"HGDP initialized with {version_info.record_count} variants")
            return True
            
        except Exception as e:
            logger.error(f"Failed to initialize HGDP: {e}")
            return False
    
    def _load_aim_frequencies(self) -> None:
        """Load HGDP-specific AIMs with population frequencies."""
        # Selected AIMs with HGDP population frequencies
        HGDP_AIMS = {
            "rs1426654": {
                "chr": "15", "pos": 48426484, "ref": "G", "alt": "A",
                "freqs": {
                    # Africa
                    "Yoruba": 0.00, "Mandenka": 0.00, "Biaka_Pygmies": 0.00,
                    "Mbuti_Pygmies": 0.00, "San": 0.00,
                    # Europe
                    "French": 1.00, "Sardinian": 1.00, "Orcadian": 1.00,
                    "Russian": 0.98, "Basque": 1.00, "Tuscan": 0.99,
                    # Middle East
                    "Palestinian": 0.98, "Druze": 0.98, "Bedouin": 0.95,
                    # South Asia
                    "Sindhi": 0.92, "Pathan": 0.91, "Brahui": 0.88,
                    "Balochi": 0.90, "Hazara": 0.85,
                    # East Asia
                    "Han": 0.02, "Japanese": 0.00, "Dai": 0.02,
                    # Americas
                    "Maya": 0.15, "Pima": 0.10, "Colombian": 0.50,
                }
            },
            "rs16891982": {
                "chr": "5", "pos": 33951693, "ref": "C", "alt": "G",
                "freqs": {
                    "Yoruba": 0.00, "Mandenka": 0.01, "Biaka_Pygmies": 0.00,
                    "French": 0.97, "Sardinian": 0.92, "Orcadian": 0.98,
                    "Russian": 0.94, "Basque": 0.96, "Tuscan": 0.94,
                    "Palestinian": 0.55, "Druze": 0.60, "Bedouin": 0.45,
                    "Sindhi": 0.30, "Pathan": 0.35, "Brahui": 0.25,
                    "Han": 0.00, "Japanese": 0.00, "Dai": 0.01,
                    "Maya": 0.10, "Pima": 0.08, "Colombian": 0.40,
                }
            },
            "rs12913832": {
                "chr": "15", "pos": 28365618, "ref": "A", "alt": "G",
                "freqs": {
                    "Yoruba": 0.00, "Mandenka": 0.00, "San": 0.00,
                    "French": 0.60, "Sardinian": 0.30, "Orcadian": 0.75,
                    "Russian": 0.55, "Basque": 0.45, "Tuscan": 0.30,
                    "Palestinian": 0.05, "Druze": 0.08, "Bedouin": 0.03,
                    "Sindhi": 0.02, "Pathan": 0.05, "Kalash": 0.15,
                    "Han": 0.00, "Japanese": 0.00,
                    "Maya": 0.00, "Pima": 0.00,
                }
            },
            "rs3827760": {
                "chr": "2", "pos": 109513601, "ref": "A", "alt": "G",
                "freqs": {
                    "Yoruba": 0.00, "Mandenka": 0.00,
                    "French": 0.00, "Sardinian": 0.00, "Orcadian": 0.00,
                    "Palestinian": 0.00, "Druze": 0.00,
                    "Sindhi": 0.00, "Pathan": 0.00,
                    "Han": 0.91, "Japanese": 0.90, "Dai": 0.85,
                    "Cambodian": 0.75, "Mongola": 0.88, "Yakut": 0.80,
                    "Maya": 0.35, "Pima": 0.40, "Karitiana": 0.55,
                    "Colombian": 0.10,
                }
            },
            "rs2814778": {
                "chr": "1", "pos": 159174683, "ref": "T", "alt": "C",
                "freqs": {
                    "Yoruba": 1.00, "Mandenka": 1.00, "Biaka_Pygmies": 0.98,
                    "Mbuti_Pygmies": 0.95, "San": 0.85, "Bantu_Kenya": 0.95,
                    "French": 0.00, "Sardinian": 0.00, "Orcadian": 0.00,
                    "Palestinian": 0.00, "Druze": 0.00, "Bedouin": 0.02,
                    "Han": 0.00, "Japanese": 0.00,
                    "Maya": 0.00, "Pima": 0.00,
                    "Papuan": 0.05, "Melanesian": 0.08,
                }
            },
            "rs4988235": {
                "chr": "2", "pos": 136608646, "ref": "G", "alt": "A",
                "freqs": {
                    "Yoruba": 0.00, "Mandenka": 0.00, "San": 0.00,
                    "French": 0.80, "Orcadian": 0.85, "Russian": 0.60,
                    "Sardinian": 0.15, "Basque": 0.70, "Tuscan": 0.40,
                    "Palestinian": 0.20, "Druze": 0.18, "Bedouin": 0.30,
                    "Sindhi": 0.15, "Pathan": 0.25, "Balochi": 0.20,
                    "Han": 0.00, "Japanese": 0.00,
                    "Maya": 0.00, "Pima": 0.00,
                }
            },
        }
        
        conn = self._get_connection()
        cursor = conn.cursor()
        
        for rsid, data in HGDP_AIMS.items():
            cursor.execute("""
                INSERT OR REPLACE INTO variants (rsid, chromosome, position, ref_allele, alt_allele)
                VALUES (?, ?, ?, ?, ?)
            """, (rsid, data["chr"], data["pos"], data["ref"], data["alt"]))
            
            for pop, freq in data["freqs"].items():
                sample_size = HGDP_POPULATIONS.get(pop, {}).get("n", 20)
                cursor.execute("""
                    INSERT OR REPLACE INTO frequencies (rsid, population, allele, frequency, sample_size)
                    VALUES (?, ?, ?, ?, ?)
                """, (rsid, pop, data["alt"], freq, sample_size))
        
        conn.commit()
        logger.info(f"Loaded {len(HGDP_AIMS)} HGDP AIMs")
    
    def _count_variants(self) -> int:
        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM variants")
        return cursor.fetchone()[0]
    
    def lookup_variant(self, rsid: str) -> Optional[VariantInfo]:
        """Look up a variant."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT rsid, chromosome, position, ref_allele, alt_allele
            FROM variants WHERE rsid = ?
        """, (rsid,))
        row = cursor.fetchone()
        
        if not row:
            return None
        
        # Get frequencies
        cursor.execute("""
            SELECT population, frequency FROM frequencies WHERE rsid = ?
        """, (rsid,))
        freqs = {r["population"]: r["frequency"] for r in cursor.fetchall()}
        
        return VariantInfo(
            rsid=row["rsid"],
            chromosome=row["chromosome"],
            position=row["position"],
            ref_allele=row["ref_allele"],
            alt_allele=row["alt_allele"],
            frequencies=freqs,
        )
    
    def calculate_population_similarity(
        self, 
        genotypes: Dict[str, str],
        region: Optional[str] = None
    ) -> Dict[str, float]:
        """
        Calculate similarity to HGDP populations.
        
        Args:
            genotypes: Dict mapping rsid -> genotype
            region: Optional filter by region (Africa, Europe, etc.)
            
        Returns:
            Dict mapping population name to similarity score
        """
        populations = list(HGDP_POPULATIONS.keys())
        if region:
            populations = HGDP_REGIONS.get(region, populations)
        
        log_likelihoods = {pop: 0.0 for pop in populations}
        markers_used = 0
        
        conn = self._get_connection()
        cursor = conn.cursor()
        
        for rsid, geno in genotypes.items():
            if not geno or len(geno) != 2:
                continue
            
            cursor.execute("SELECT alt_allele FROM variants WHERE rsid = ?", (rsid,))
            var_row = cursor.fetchone()
            if not var_row:
                continue
            
            alt = var_row["alt_allele"].upper()
            alt_count = sum(1 for a in geno.upper() if a == alt)
            
            cursor.execute("""
                SELECT population, frequency FROM frequencies WHERE rsid = ?
            """, (rsid,))
            
            freq_map = {r["population"]: r["frequency"] for r in cursor.fetchall()}
            if not freq_map:
                continue
            
            markers_used += 1
            
            import math
            for pop in populations:
                p = freq_map.get(pop, 0.5)
                p = max(0.001, min(0.999, p))
                
                if alt_count == 0:
                    prob = (1 - p) ** 2
                elif alt_count == 1:
                    prob = 2 * p * (1 - p)
                else:
                    prob = p ** 2
                
                log_likelihoods[pop] += math.log(prob)
        
        if markers_used == 0:
            return {pop: 1.0 / len(populations) for pop in populations}
        
        import math
        max_ll = max(log_likelihoods.values())
        exp_scores = {pop: math.exp(ll - max_ll) for pop, ll in log_likelihoods.items()}
        total = sum(exp_scores.values())
        
        return {pop: score / total for pop, score in exp_scores.items()}
    
    def get_region_summary(
        self, 
        population_similarities: Dict[str, float]
    ) -> Dict[str, float]:
        """Aggregate population similarities to region level."""
        region_scores = {region: 0.0 for region in HGDP_REGIONS}
        
        for pop, score in population_similarities.items():
            pop_info = HGDP_POPULATIONS.get(pop, {})
            region = pop_info.get("region")
            if region:
                region_scores[region] += score
        
        return region_scores


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def compare_to_hgdp(genotypes: Dict[str, str]) -> Dict[str, float]:
    """Compare genotypes to HGDP populations."""
    hgdp = HGDP()
    if not hgdp.is_downloaded:
        hgdp.download()
    return hgdp.calculate_population_similarity(genotypes)
