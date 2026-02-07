"""
Base classes and utilities for genomics reference datasets.

All datasets are downloaded and cached locally for privacy-first analysis.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union
import hashlib
import json
import gzip
import logging
import os
import urllib.request
import ssl

logger = logging.getLogger(__name__)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Base path for all reference datasets
DATASETS_BASE_PATH = Path.home() / ".openclaw" / "workspace" / "skills" / "personal-genomics" / "references" / "datasets"


@dataclass
class DatasetVersion:
    """Track dataset version information."""
    name: str
    version: str
    downloaded: datetime
    source_url: str
    file_hash: Optional[str] = None
    record_count: Optional[int] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "version": self.version,
            "downloaded": self.downloaded.isoformat(),
            "source_url": self.source_url,
            "file_hash": self.file_hash,
            "record_count": self.record_count,
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DatasetVersion":
        return cls(
            name=data["name"],
            version=data["version"],
            downloaded=datetime.fromisoformat(data["downloaded"]),
            source_url=data["source_url"],
            file_hash=data.get("file_hash"),
            record_count=data.get("record_count"),
        )


@dataclass
class VariantInfo:
    """Standard variant information structure."""
    rsid: str
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    gene: Optional[str] = None
    consequence: Optional[str] = None
    
    # Population frequencies
    frequencies: Dict[str, float] = field(default_factory=dict)
    
    # Clinical info (if applicable)
    clinical_significance: Optional[str] = None
    conditions: List[str] = field(default_factory=list)
    
    # References
    pmids: List[str] = field(default_factory=list)


@dataclass 
class PopulationFrequency:
    """Allele frequency in a reference population."""
    population: str
    population_name: str
    superpopulation: Optional[str]
    allele: str
    frequency: float
    allele_count: int
    total_alleles: int
    
    @property
    def sample_size(self) -> int:
        return self.total_alleles // 2


# =============================================================================
# BASE DATASET CLASS
# =============================================================================

class BaseDataset(ABC):
    """
    Base class for all reference datasets.
    
    Provides common functionality for:
    - Downloading and caching data
    - Version tracking
    - Lookup methods
    """
    
    name: str = "base"
    version: str = "unknown"
    description: str = ""
    source_url: str = ""
    
    def __init__(self, data_dir: Optional[Path] = None):
        self.data_dir = data_dir or DATASETS_BASE_PATH / self.name
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self._version_info: Optional[DatasetVersion] = None
        self._cache: Dict[str, Any] = {}
        
    @property
    def version_file(self) -> Path:
        return self.data_dir / "version.json"
    
    @property
    def is_downloaded(self) -> bool:
        """Check if dataset has been downloaded."""
        return self.version_file.exists()
    
    def get_version_info(self) -> Optional[DatasetVersion]:
        """Get current version information."""
        if self._version_info is None and self.version_file.exists():
            with open(self.version_file) as f:
                data = json.load(f)
                self._version_info = DatasetVersion.from_dict(data)
        return self._version_info
    
    def save_version_info(self, version_info: DatasetVersion) -> None:
        """Save version information."""
        with open(self.version_file, "w") as f:
            json.dump(version_info.to_dict(), f, indent=2)
        self._version_info = version_info
    
    @abstractmethod
    def download(self, force: bool = False) -> bool:
        """
        Download the dataset.
        
        Args:
            force: Re-download even if already present
            
        Returns:
            True if download successful
        """
        pass
    
    @abstractmethod
    def lookup_variant(self, rsid: str) -> Optional[VariantInfo]:
        """
        Look up a variant by rsID.
        
        Args:
            rsid: The variant identifier (e.g., rs1426654)
            
        Returns:
            VariantInfo if found, None otherwise
        """
        pass
    
    def lookup_variants(self, rsids: List[str]) -> Dict[str, Optional[VariantInfo]]:
        """
        Look up multiple variants.
        
        Args:
            rsids: List of variant identifiers
            
        Returns:
            Dict mapping rsid to VariantInfo (or None if not found)
        """
        return {rsid: self.lookup_variant(rsid) for rsid in rsids}
    
    def _download_file(
        self, 
        url: str, 
        dest: Path, 
        description: str = "file"
    ) -> bool:
        """
        Download a file with progress logging.
        
        Args:
            url: URL to download
            dest: Destination path
            description: Description for logging
            
        Returns:
            True if successful
        """
        logger.info(f"Downloading {description} from {url}")
        
        # Create SSL context that doesn't verify (some FTP servers have cert issues)
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        
        try:
            # Handle both HTTP and FTP
            if url.startswith("ftp://"):
                urllib.request.urlretrieve(url, dest)
            else:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, context=ctx, timeout=300) as response:
                    with open(dest, "wb") as out:
                        out.write(response.read())
            
            logger.info(f"Downloaded {description} to {dest}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to download {description}: {e}")
            if dest.exists():
                dest.unlink()
            return False
    
    def _compute_file_hash(self, filepath: Path) -> str:
        """Compute MD5 hash of a file."""
        hash_md5 = hashlib.md5()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()


# =============================================================================
# INDEX-BASED DATASET (for large datasets using rsID index files)
# =============================================================================

class IndexedDataset(BaseDataset):
    """
    Base class for large datasets using index files for fast lookups.
    
    Uses a two-level index:
    1. rsID -> file offset (stored in index file)
    2. Seek to offset to read record
    """
    
    @property
    def index_file(self) -> Path:
        return self.data_dir / "rsid_index.json"
    
    @property
    def data_file(self) -> Path:
        return self.data_dir / "data.tsv.gz"
    
    def _build_index(self, data_file: Path, rsid_col: int = 0) -> Dict[str, int]:
        """
        Build an index mapping rsID to line number.
        
        For compressed files, we index by line number since seeking
        in gzip files is not efficient.
        """
        index = {}
        line_num = 0
        
        opener = gzip.open if str(data_file).endswith(".gz") else open
        mode = "rt" if str(data_file).endswith(".gz") else "r"
        
        with opener(data_file, mode) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) > rsid_col:
                    rsid = parts[rsid_col]
                    if rsid.startswith("rs"):
                        index[rsid] = line_num
                line_num += 1
        
        return index
    
    def _load_index(self) -> Dict[str, int]:
        """Load the rsID index from disk."""
        if "index" not in self._cache:
            if self.index_file.exists():
                with open(self.index_file) as f:
                    self._cache["index"] = json.load(f)
            else:
                self._cache["index"] = {}
        return self._cache["index"]
    
    def _save_index(self, index: Dict[str, int]) -> None:
        """Save the rsID index to disk."""
        with open(self.index_file, "w") as f:
            json.dump(index, f)
        self._cache["index"] = index


# =============================================================================
# SQLITE-BASED DATASET (for complex queries)
# =============================================================================

class SQLiteDataset(BaseDataset):
    """
    Base class for datasets using SQLite for storage.
    
    Better for complex queries and large datasets with multiple lookup patterns.
    """
    
    @property
    def db_file(self) -> Path:
        return self.data_dir / "data.db"
    
    def _get_connection(self):
        """Get SQLite connection."""
        import sqlite3
        if "conn" not in self._cache:
            self._cache["conn"] = sqlite3.connect(str(self.db_file))
            self._cache["conn"].row_factory = sqlite3.Row
        return self._cache["conn"]
    
    def close(self):
        """Close database connection."""
        if "conn" in self._cache:
            self._cache["conn"].close()
            del self._cache["conn"]


# =============================================================================
# POPULATION CODES
# =============================================================================

# 1000 Genomes population codes
THOUSAND_GENOMES_POPULATIONS = {
    # European (EUR)
    "CEU": {"name": "Utah residents with Northern/Western European ancestry", "superpop": "EUR"},
    "GBR": {"name": "British in England and Scotland", "superpop": "EUR"},
    "FIN": {"name": "Finnish in Finland", "superpop": "EUR"},
    "IBS": {"name": "Iberian Population in Spain", "superpop": "EUR"},
    "TSI": {"name": "Toscani in Italia", "superpop": "EUR"},
    
    # African (AFR)
    "YRI": {"name": "Yoruba in Ibadan, Nigeria", "superpop": "AFR"},
    "LWK": {"name": "Luhya in Webuye, Kenya", "superpop": "AFR"},
    "GWD": {"name": "Gambian in Western Divisions in the Gambia", "superpop": "AFR"},
    "MSL": {"name": "Mende in Sierra Leone", "superpop": "AFR"},
    "ESN": {"name": "Esan in Nigeria", "superpop": "AFR"},
    "ASW": {"name": "African Ancestry in Southwest US", "superpop": "AFR"},
    "ACB": {"name": "African Caribbean in Barbados", "superpop": "AFR"},
    
    # East Asian (EAS)
    "CHB": {"name": "Han Chinese in Beijing, China", "superpop": "EAS"},
    "JPT": {"name": "Japanese in Tokyo, Japan", "superpop": "EAS"},
    "CHS": {"name": "Southern Han Chinese", "superpop": "EAS"},
    "CDX": {"name": "Chinese Dai in Xishuangbanna, China", "superpop": "EAS"},
    "KHV": {"name": "Kinh in Ho Chi Minh City, Vietnam", "superpop": "EAS"},
    
    # South Asian (SAS)
    "GIH": {"name": "Gujarati Indian in Houston, Texas", "superpop": "SAS"},
    "PJL": {"name": "Punjabi in Lahore, Pakistan", "superpop": "SAS"},
    "BEB": {"name": "Bengali in Bangladesh", "superpop": "SAS"},
    "STU": {"name": "Sri Lankan Tamil in the UK", "superpop": "SAS"},
    "ITU": {"name": "Indian Telugu in the UK", "superpop": "SAS"},
    
    # Americas (AMR)
    "MXL": {"name": "Mexican Ancestry in Los Angeles", "superpop": "AMR"},
    "PUR": {"name": "Puerto Rican in Puerto Rico", "superpop": "AMR"},
    "CLM": {"name": "Colombian in Medellin, Colombia", "superpop": "AMR"},
    "PEL": {"name": "Peruvian in Lima, Peru", "superpop": "AMR"},
}

SUPERPOPULATIONS = {
    "EUR": "European",
    "AFR": "African",
    "EAS": "East Asian",
    "SAS": "South Asian",
    "AMR": "Ad Mixed American",
}


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_all_datasets() -> Dict[str, type]:
    """Get all available dataset classes."""
    from . import (
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
    return {
        "1000genomes": ThousandGenomes,
        "gnomad": GnomAD,
        "clinvar": ClinVar,
        "pharmgkb": PharmGKB,
        "pgs_catalog": PGSCatalog,
        "gwas_catalog": GWASCatalog,
        "hgdp": HGDP,
        "sgdp": SGDP,
        "ancient_dna": AncientDNA,
    }


def download_all_datasets(force: bool = False) -> Dict[str, bool]:
    """Download all datasets."""
    results = {}
    for name, cls in get_all_datasets().items():
        try:
            dataset = cls()
            results[name] = dataset.download(force=force)
        except Exception as e:
            logger.error(f"Failed to download {name}: {e}")
            results[name] = False
    return results


def get_dataset_status() -> Dict[str, Dict[str, Any]]:
    """Get status of all datasets."""
    status = {}
    for name, cls in get_all_datasets().items():
        try:
            dataset = cls()
            version_info = dataset.get_version_info()
            status[name] = {
                "downloaded": dataset.is_downloaded,
                "version": version_info.version if version_info else None,
                "downloaded_date": version_info.downloaded.isoformat() if version_info else None,
                "record_count": version_info.record_count if version_info else None,
            }
        except Exception as e:
            status[name] = {"error": str(e)}
    return status
