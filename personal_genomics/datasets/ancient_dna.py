"""
Ancient DNA Dataset Integration

Provides markers for detecting ancestral population signals:
- Western Hunter-Gatherers (WHG)
- Eastern Hunter-Gatherers (EHG)
- Scandinavian Hunter-Gatherers (SHG)
- Neolithic Farmers (Anatolian/European)
- Yamnaya/Steppe Pastoralists
- Neanderthal admixture
- Denisovan admixture

Data sources:
- Reich Lab public ancient DNA (https://reich.hms.harvard.edu/datasets)
- Haak et al. 2015 - Massive migration from the steppe
- Mathieson et al. 2015 - Genome-wide patterns of selection
- Lazaridis et al. 2014 - Ancient human genomes
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import json
import logging

from .base import (
    BaseDataset,
    DatasetVersion,
    VariantInfo,
    DATASETS_BASE_PATH,
)

logger = logging.getLogger(__name__)


# =============================================================================
# ANCIENT POPULATION MARKERS
# =============================================================================

@dataclass
class AncientMarker:
    """A marker informative for ancient population ancestry."""
    rsid: str
    chromosome: str
    position: int
    ancestral_allele: str
    derived_allele: str
    population: str  # WHG, EHG, Neolithic, Yamnaya, Neanderthal, Denisovan
    derived_frequency: float  # Frequency of derived allele in ancient population
    phenotype: Optional[str] = None  # Associated phenotype if known
    pmid: Optional[str] = None  # PubMed reference
    notes: Optional[str] = None


# Curated ancestry-informative markers from published ancient DNA studies
# These are well-established markers with clear population-specific signals
ANCIENT_MARKERS: List[AncientMarker] = [
    # =========================================================================
    # WESTERN HUNTER-GATHERER (WHG) MARKERS
    # Pre-Neolithic Europeans, ~14,000-5,000 BCE
    # =========================================================================
    AncientMarker(
        rsid="rs12913832",
        chromosome="15",
        position=28365618,
        ancestral_allele="A",
        derived_allele="G",
        population="WHG",
        derived_frequency=0.95,
        phenotype="Blue eyes (HERC2/OCA2)",
        pmid="18172690",
        notes="Nearly fixed in WHG, signature of European hunter-gatherers"
    ),
    AncientMarker(
        rsid="rs1800407",
        chromosome="15",
        position=28230318,
        ancestral_allele="G",
        derived_allele="A",
        population="WHG",
        derived_frequency=0.85,
        phenotype="Blue/green eyes (OCA2)",
        pmid="17952075",
        notes="High frequency in Mesolithic Europeans"
    ),
    AncientMarker(
        rsid="rs16891982",
        chromosome="5",
        position=33951693,
        ancestral_allele="C",
        derived_allele="G",
        population="WHG",
        derived_frequency=0.80,
        phenotype="Light skin (SLC45A2)",
        pmid="17182896",
        notes="Present in WHG but not fixed like in later Europeans"
    ),
    AncientMarker(
        rsid="rs1426654",
        chromosome="15",
        position=48426484,
        ancestral_allele="A",
        derived_allele="G",
        population="WHG",
        derived_frequency=0.20,
        phenotype="Light skin (SLC24A5)",
        pmid="16357253",
        notes="LOW in WHG - this came with Neolithic farmers"
    ),

    # =========================================================================
    # NEOLITHIC FARMER MARKERS
    # Anatolian/Early European Farmers, ~9,000-4,000 BCE
    # =========================================================================
    AncientMarker(
        rsid="rs1426654",
        chromosome="15",
        position=48426484,
        ancestral_allele="A",
        derived_allele="G",
        population="Neolithic_Farmer",
        derived_frequency=0.95,
        phenotype="Light skin (SLC24A5)",
        pmid="16357253",
        notes="Nearly fixed in Neolithic farmers, brought to Europe"
    ),
    AncientMarker(
        rsid="rs4988235",
        chromosome="2",
        position=136608646,
        ancestral_allele="G",
        derived_allele="A",
        population="Neolithic_Farmer",
        derived_frequency=0.10,
        phenotype="Lactase persistence (LCT)",
        pmid="11788828",
        notes="RARE in Neolithic farmers - selected later in Bronze Age"
    ),
    AncientMarker(
        rsid="rs2282679",
        chromosome="4",
        position=72618323,
        ancestral_allele="T",
        derived_allele="G",
        population="Neolithic_Farmer",
        derived_frequency=0.65,
        phenotype="Lower vitamin D binding protein (GC)",
        pmid="20541252",
        notes="Associated with vitamin D levels, higher in farmers"
    ),

    # =========================================================================
    # YAMNAYA / STEPPE PASTORALIST MARKERS
    # Pontic-Caspian Steppe, ~3,300-2,600 BCE
    # =========================================================================
    AncientMarker(
        rsid="rs4988235",
        chromosome="2",
        position=136608646,
        ancestral_allele="G",
        derived_allele="A",
        population="Yamnaya",
        derived_frequency=0.25,
        phenotype="Lactase persistence (LCT)",
        pmid="11788828",
        notes="Present in Yamnaya, underwent strong selection"
    ),
    AncientMarker(
        rsid="rs16891982",
        chromosome="5",
        position=33951693,
        ancestral_allele="C",
        derived_allele="G",
        population="Yamnaya",
        derived_frequency=0.90,
        phenotype="Light skin (SLC45A2)",
        pmid="17182896",
        notes="High frequency in Steppe populations"
    ),
    AncientMarker(
        rsid="rs12821256",
        chromosome="12",
        position=89328335,
        ancestral_allele="T",
        derived_allele="C",
        population="Yamnaya",
        derived_frequency=0.30,
        phenotype="Blond hair (KITLG)",
        pmid="24531881",
        notes="Associated with lighter hair, moderate in Steppe"
    ),
    AncientMarker(
        rsid="rs2814778",
        chromosome="1",
        position=159174683,
        ancestral_allele="T",
        derived_allele="C",
        population="Yamnaya",
        derived_frequency=0.98,
        phenotype="Duffy-null (DARC)",
        pmid="11509864",
        notes="Fixed ancestral in all non-African populations"
    ),

    # =========================================================================
    # EASTERN HUNTER-GATHERER (EHG) MARKERS
    # Eastern European / Western Siberian, ~10,000-5,000 BCE
    # =========================================================================
    AncientMarker(
        rsid="rs12913832",
        chromosome="15",
        position=28365618,
        ancestral_allele="A",
        derived_allele="G",
        population="EHG",
        derived_frequency=0.70,
        phenotype="Blue eyes (HERC2/OCA2)",
        pmid="18172690",
        notes="Present but not as fixed as in WHG"
    ),
    AncientMarker(
        rsid="rs1042602",
        chromosome="11",
        position=88911696,
        ancestral_allele="C",
        derived_allele="A",
        population="EHG",
        derived_frequency=0.60,
        phenotype="Lighter skin (TYR)",
        pmid="20421989",
        notes="Moderate frequency in EHG"
    ),

    # =========================================================================
    # NEANDERTHAL INTROGRESSION MARKERS
    # Archaic admixture, ~50,000-40,000 BCE
    # =========================================================================
    AncientMarker(
        rsid="rs2298813",
        chromosome="16",
        position=89985940,
        ancestral_allele="A",
        derived_allele="G",
        population="Neanderthal",
        derived_frequency=0.70,
        phenotype="BNC2 - skin pigmentation",
        pmid="26140447",
        notes="Neanderthal introgressed allele affecting skin"
    ),
    AncientMarker(
        rsid="rs10519177",
        chromosome="9",
        position=16802314,
        ancestral_allele="G",
        derived_allele="A",
        population="Neanderthal",
        derived_frequency=0.65,
        phenotype="BNC2 region",
        pmid="26140447",
        notes="Linked to Neanderthal introgression"
    ),
    AncientMarker(
        rsid="rs2066853",
        chromosome="1",
        position=17345517,
        ancestral_allele="G",
        derived_allele="A",
        population="Neanderthal",
        derived_frequency=0.50,
        phenotype="AHR - immune response",
        pmid="26417903",
        notes="Neanderthal variant affecting immune function"
    ),
    AncientMarker(
        rsid="rs11803731",
        chromosome="1",
        position=21833535,
        ancestral_allele="A",
        derived_allele="T",
        population="Neanderthal",
        derived_frequency=0.55,
        phenotype="TCHH - hair morphology",
        pmid="26417903",
        notes="Neanderthal introgressed, affects hair texture"
    ),

    # =========================================================================
    # DENISOVAN INTROGRESSION MARKERS
    # Primarily detected in East Asian/Oceanian populations
    # =========================================================================
    AncientMarker(
        rsid="rs1058172",
        chromosome="2",
        position=109513601,
        ancestral_allele="C",
        derived_allele="T",
        population="Denisovan",
        derived_frequency=0.40,
        phenotype="EDAR - hair/teeth/sweat glands",
        pmid="23151665",
        notes="East Asian selection, possibly Denisovan-linked"
    ),
    AncientMarker(
        rsid="rs3827760",
        chromosome="2",
        position=109513601,
        ancestral_allele="A",
        derived_allele="G",
        population="Denisovan",
        derived_frequency=0.85,
        phenotype="EDAR V370A - thick hair, shovel incisors",
        pmid="23151665",
        notes="Strong East Asian signal, archaic introgression"
    ),
]

# Group markers by population for easy lookup
MARKERS_BY_POPULATION: Dict[str, List[AncientMarker]] = {}
for marker in ANCIENT_MARKERS:
    if marker.population not in MARKERS_BY_POPULATION:
        MARKERS_BY_POPULATION[marker.population] = []
    MARKERS_BY_POPULATION[marker.population].append(marker)


@dataclass
class AncestralSignal:
    """Result of ancestral population signal detection."""
    population: str
    description: str
    markers_tested: int
    markers_found: int
    derived_alleles: int
    signal_strength: str  # "Strong", "Moderate", "Weak", "None"
    confidence: str  # "High", "Medium", "Low"
    phenotypic_evidence: List[str]
    notes: str


class AncientDNADataset(BaseDataset):
    """
    Ancient DNA reference dataset for detecting ancestral population signals.
    
    Uses curated markers from published ancient DNA studies to detect:
    - Western/Eastern Hunter-Gatherer ancestry
    - Neolithic Farmer ancestry  
    - Yamnaya/Steppe ancestry
    - Neanderthal/Denisovan introgression
    """
    
    name = "ancient_dna"
    version = "1.0.0"
    description = "Ancient DNA ancestral markers"
    source_url = "curated_from_literature"
    
    def __init__(self, data_dir: Optional[Path] = None):
        super().__init__(data_dir or DATASETS_BASE_PATH / "ancient_dna")
        self._markers = {m.rsid: m for m in ANCIENT_MARKERS}
        
    def download(self, force: bool = False) -> bool:
        """
        No download needed - markers are curated from published literature.
        This method exists for API consistency.
        """
        if self.is_downloaded and not force:
            logger.info("Ancient DNA markers already loaded")
            return True
            
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # Save marker database to JSON for transparency
        markers_file = self.data_dir / "curated_markers.json"
        markers_data = {
            "version": self.version,
            "source": "Curated from published ancient DNA studies",
            "references": [
                "Haak et al. 2015 (PMID: 25731166)",
                "Mathieson et al. 2015 (PMID: 26595274)",
                "Lazaridis et al. 2014 (PMID: 25230663)",
                "Sankararaman et al. 2014 (PMID: 24476815)",
            ],
            "marker_count": len(ANCIENT_MARKERS),
            "populations": list(MARKERS_BY_POPULATION.keys()),
            "markers": [
                {
                    "rsid": m.rsid,
                    "chromosome": m.chromosome,
                    "position": m.position,
                    "population": m.population,
                    "phenotype": m.phenotype,
                    "pmid": m.pmid,
                }
                for m in ANCIENT_MARKERS
            ]
        }
        
        with open(markers_file, "w") as f:
            json.dump(markers_data, f, indent=2)
            
        # Save version info
        version_info = DatasetVersion(
            name=self.name,
            version=self.version,
            downloaded=datetime.now(),
            source_url=self.source_url,
            record_count=len(ANCIENT_MARKERS),
        )
        self.save_version_info(version_info)
        
        return True
    
    def is_available(self) -> bool:
        """Always available - markers are embedded."""
        return True
    
    def lookup_variant(self, rsid: str) -> Optional[VariantInfo]:
        """Look up a variant by rsID (required by BaseDataset)."""
        marker = self._markers.get(rsid)
        if marker:
            return VariantInfo(
                rsid=marker.rsid,
                chromosome=marker.chromosome,
                position=marker.position,
                ref_allele=marker.ancestral_allele,
                alt_allele=marker.derived_allele,
            )
        return None
    
    def get_marker(self, rsid: str) -> Optional[AncientMarker]:
        """Get marker info by rsID."""
        return self._markers.get(rsid)
    
    def get_population_markers(self, population: str) -> List[AncientMarker]:
        """Get all markers for a specific ancient population."""
        return MARKERS_BY_POPULATION.get(population, [])
    
    def get_all_populations(self) -> List[str]:
        """Get list of all ancient populations."""
        return list(MARKERS_BY_POPULATION.keys())
    
    def analyze_ancestral_signals(
        self, 
        genotypes: Dict[str, str]
    ) -> Dict[str, AncestralSignal]:
        """
        Analyze genotypes for ancestral population signals.
        
        Args:
            genotypes: Dict mapping rsID -> genotype (e.g., "rs12913832" -> "AG")
            
        Returns:
            Dict mapping population name -> AncestralSignal result
        """
        results = {}
        
        for population, markers in MARKERS_BY_POPULATION.items():
            markers_tested = 0
            markers_found = 0
            derived_count = 0
            phenotypes_found = []
            
            for marker in markers:
                if marker.rsid in genotypes:
                    markers_tested += 1
                    markers_found += 1
                    genotype = genotypes[marker.rsid].upper()
                    
                    # Count derived alleles
                    derived = marker.derived_allele.upper()
                    derived_in_geno = genotype.count(derived)
                    derived_count += derived_in_geno
                    
                    if derived_in_geno > 0 and marker.phenotype:
                        phenotypes_found.append(marker.phenotype)
            
            # Calculate signal strength
            if markers_found == 0:
                signal_strength = "None"
                confidence = "Low"
                notes = f"No {population} markers found in data"
            else:
                # Normalize: max possible derived is 2 per marker
                max_derived = markers_found * 2
                derived_ratio = derived_count / max_derived if max_derived > 0 else 0
                
                if derived_ratio >= 0.7:
                    signal_strength = "Strong"
                    confidence = "High" if markers_found >= 3 else "Medium"
                elif derived_ratio >= 0.4:
                    signal_strength = "Moderate"
                    confidence = "Medium" if markers_found >= 2 else "Low"
                elif derived_ratio >= 0.2:
                    signal_strength = "Weak"
                    confidence = "Low"
                else:
                    signal_strength = "None"
                    confidence = "Low"
                
                notes = f"{derived_count}/{max_derived} derived alleles across {markers_found} markers"
            
            # Population descriptions
            descriptions = {
                "WHG": "Western Hunter-Gatherers: Pre-Neolithic Europeans (~14,000-5,000 BCE)",
                "EHG": "Eastern Hunter-Gatherers: Eastern European/Western Siberian (~10,000-5,000 BCE)",
                "Neolithic_Farmer": "Neolithic Farmers: Anatolian/Early European Farmers (~9,000-4,000 BCE)",
                "Yamnaya": "Yamnaya/Steppe Pastoralists: Pontic-Caspian Steppe (~3,300-2,600 BCE)",
                "Neanderthal": "Neanderthal Introgression: Archaic human admixture (~50,000-40,000 BCE)",
                "Denisovan": "Denisovan Introgression: Archaic human admixture (primarily East Asian/Oceanian)",
            }
            
            results[population] = AncestralSignal(
                population=population,
                description=descriptions.get(population, population),
                markers_tested=markers_tested,
                markers_found=markers_found,
                derived_alleles=derived_count,
                signal_strength=signal_strength,
                confidence=confidence,
                phenotypic_evidence=list(set(phenotypes_found)),
                notes=notes,
            )
        
        return results


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def get_ancient_markers() -> List[AncientMarker]:
    """Get all curated ancient DNA markers."""
    return ANCIENT_MARKERS.copy()


def analyze_ancient_ancestry(genotypes: Dict[str, str]) -> Dict[str, AncestralSignal]:
    """
    Convenience function to analyze ancestral signals.
    
    Args:
        genotypes: Dict mapping rsID -> genotype
        
    Returns:
        Dict of population -> AncestralSignal
    """
    dataset = AncientDNADataset()
    return dataset.analyze_ancestral_signals(genotypes)


def get_european_ancestry_narrative(signals: Dict[str, AncestralSignal]) -> str:
    """
    Generate a narrative description of European ancestry from ancient signals.
    
    Most modern Europeans are a mix of:
    1. Western Hunter-Gatherers (WHG) - ~10-20%
    2. Neolithic Farmers (Anatolian) - ~30-50%  
    3. Yamnaya/Steppe - ~30-50%
    
    Plus ~2% Neanderthal admixture.
    """
    narrative_parts = []
    
    # Check for WHG signal
    if "WHG" in signals and signals["WHG"].signal_strength in ["Strong", "Moderate"]:
        narrative_parts.append(
            f"**Western Hunter-Gatherer Signal ({signals['WHG'].signal_strength})**: "
            f"You carry markers associated with Mesolithic Europeans who lived in Western Europe "
            f"before the arrival of farming (~14,000-5,000 BCE). "
            f"Phenotypic evidence: {', '.join(signals['WHG'].phenotypic_evidence) or 'None detected'}."
        )
    
    # Check for Neolithic Farmer signal
    if "Neolithic_Farmer" in signals and signals["Neolithic_Farmer"].signal_strength in ["Strong", "Moderate"]:
        narrative_parts.append(
            f"**Neolithic Farmer Signal ({signals['Neolithic_Farmer'].signal_strength})**: "
            f"You carry markers from the Neolithic farmers who migrated from Anatolia/Near East "
            f"and brought agriculture to Europe (~9,000-4,000 BCE). "
            f"Phenotypic evidence: {', '.join(signals['Neolithic_Farmer'].phenotypic_evidence) or 'None detected'}."
        )
    
    # Check for Yamnaya/Steppe signal
    if "Yamnaya" in signals and signals["Yamnaya"].signal_strength in ["Strong", "Moderate"]:
        narrative_parts.append(
            f"**Yamnaya/Steppe Signal ({signals['Yamnaya'].signal_strength})**: "
            f"You carry markers from the Yamnaya pastoralists who migrated from the Pontic-Caspian steppe "
            f"into Europe during the Bronze Age (~3,300-2,600 BCE). "
            f"Phenotypic evidence: {', '.join(signals['Yamnaya'].phenotypic_evidence) or 'None detected'}."
        )
    
    # Check for Neanderthal signal
    if "Neanderthal" in signals and signals["Neanderthal"].signal_strength in ["Strong", "Moderate", "Weak"]:
        narrative_parts.append(
            f"**Neanderthal Introgression ({signals['Neanderthal'].signal_strength})**: "
            f"You carry markers from Neanderthal admixture that occurred ~50,000-40,000 years ago. "
            f"Most non-African humans carry 1-4% Neanderthal DNA. "
            f"Phenotypic evidence: {', '.join(signals['Neanderthal'].phenotypic_evidence) or 'None detected'}."
        )
    
    if not narrative_parts:
        return "No strong ancient population signals detected with available markers."
    
    return "\n\n".join(narrative_parts)
