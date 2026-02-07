#!/usr/bin/env python3
"""
Ancient DNA Matching Engine
==========================

Compares user DNA to ancient individuals and cultures using Identity-by-State (IBS)
similarity. Free, open source alternative to YourTrueAncestry.

Unlike proprietary services:
- Full transparency — methodology is open and documented
- Citations for everything — every ancient sample has PMID
- Downloadable data — user can export raw comparisons
- Educational — explains what each culture/period means
- Free forever — no paywall
- Open source — code is auditable

Author: OpenClaw Personal Genomics
Version: 1.0.0
"""

import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, asdict
from collections import defaultdict
import math

# Paths
REFERENCE_DIR = Path(__file__).parent.parent / "references"
ANCIENT_INDIVIDUALS_FILE = REFERENCE_DIR / "ancient_individuals.json"
ANCIENT_CULTURES_FILE = REFERENCE_DIR / "ancient_cultures.json"


@dataclass
class AncientMatch:
    """Represents a match to an ancient individual."""
    id: str
    name: str
    similarity: float  # 0-100% (100 = identical)
    distance: float    # 0-1 (0 = identical)
    shared_snps: int
    total_snps: int
    period: str
    culture: str
    site: str
    country: str
    age_text: str
    description: str
    lat: Optional[float] = None
    lon: Optional[float] = None
    mt_haplogroup: Optional[str] = None
    y_haplogroup: Optional[str] = None
    shared_traits: Optional[List[str]] = None
    paper: Optional[str] = None
    pmid: Optional[str] = None
    
    def to_dict(self) -> Dict:
        return asdict(self)


@dataclass 
class CultureAffinity:
    """Represents affinity to an ancient culture."""
    id: str
    name: str
    affinity: float  # 0-100%
    period: str
    regions: List[str]
    description: str
    sample_count: int
    avg_shared_snps: float
    diagnostic_matches: Dict[str, str]
    map_color: str


@dataclass
class TraitComparison:
    """Comparison of specific traits between user and ancient."""
    trait_name: str
    rs_id: str
    user_genotype: str
    ancient_genotype: str
    match: bool
    trait_description: str


class AncientDNAMatcher:
    """
    Main class for ancient DNA matching.
    
    Methodology:
    -----------
    We use Identity-by-State (IBS) similarity, which counts shared alleles
    between the user and ancient samples at common SNP positions.
    
    For each SNP position:
    - If both have same homozygous genotype (e.g., AA vs AA): IBS = 2
    - If one hetero matches one allele (e.g., AG vs AA): IBS = 1
    - If both hetero (e.g., AG vs AG): IBS = 2
    - If no alleles match (e.g., AA vs GG): IBS = 0
    
    Distance = 1 - (sum of IBS scores) / (2 * number of shared SNPs)
    Similarity = 100 * (1 - distance)
    """
    
    def __init__(self):
        self.individuals = []
        self.cultures = {}
        self._load_reference_data()
    
    def _load_reference_data(self):
        """Load ancient individual and culture reference data."""
        # Load individuals
        if ANCIENT_INDIVIDUALS_FILE.exists():
            with open(ANCIENT_INDIVIDUALS_FILE) as f:
                data = json.load(f)
                self.individuals = data.get("individuals", [])
        
        # Load cultures
        if ANCIENT_CULTURES_FILE.exists():
            with open(ANCIENT_CULTURES_FILE) as f:
                data = json.load(f)
                self.cultures = data.get("cultures", {})
    
    def calculate_ibs_distance(
        self, 
        user_genotypes: Dict[str, str], 
        ancient_genotypes: Dict[str, str],
        min_snps: int = 5
    ) -> Tuple[Optional[float], int, int]:
        """
        Calculate Identity-by-State distance between user and ancient.
        
        Args:
            user_genotypes: Dict of rsid -> genotype (e.g., "rs123": "AG")
            ancient_genotypes: Dict of rsid -> genotype
            min_snps: Minimum overlapping SNPs required
            
        Returns:
            Tuple of (distance, shared_snp_count, total_ibs_score)
            distance is None if insufficient overlap
        """
        # Find shared SNPs
        shared_snps = set(user_genotypes.keys()) & set(ancient_genotypes.keys())
        
        if len(shared_snps) < min_snps:
            return (None, len(shared_snps), 0)
        
        total_ibs = 0
        valid_comparisons = 0
        
        for rs in shared_snps:
            user_geno = user_genotypes[rs].upper()
            ancient_geno = ancient_genotypes[rs].upper()
            
            # Skip if either is missing/invalid
            if len(user_geno) < 2 or len(ancient_geno) < 2:
                continue
            if user_geno in ['--', 'II', 'DD', 'DI', 'ID', '00', 'NC']:
                continue
            if ancient_geno in ['--', 'II', 'DD', 'DI', 'ID', '00', 'NC']:
                continue
            
            # Calculate IBS
            user_alleles = set(user_geno)
            ancient_alleles = set(ancient_geno)
            
            # Count matching alleles
            ibs = len(user_alleles & ancient_alleles)
            
            # If both homozygous, max IBS is 2 (both alleles match)
            # If one hetero one homo, max IBS is 1
            # If both hetero with same alleles, IBS is 2
            
            # Normalize: if genotypes identical, IBS=2
            if user_geno == ancient_geno or (set(user_geno) == set(ancient_geno)):
                ibs = 2
            elif len(user_alleles & ancient_alleles) > 0:
                # Partial match
                # AA vs AG = 1, AG vs AG = 2, AA vs GG = 0
                if len(user_alleles) == 1 and len(ancient_alleles) == 1:
                    # Both homozygous
                    ibs = 2 if user_alleles == ancient_alleles else 0
                else:
                    # At least one heterozygous
                    ibs = len(user_alleles & ancient_alleles)
            else:
                ibs = 0
            
            total_ibs += ibs
            valid_comparisons += 1
        
        if valid_comparisons < min_snps:
            return (None, valid_comparisons, total_ibs)
        
        # Distance: 0 = identical (all IBS=2), 1 = no match (all IBS=0)
        max_ibs = 2 * valid_comparisons
        distance = 1 - (total_ibs / max_ibs)
        
        return (distance, valid_comparisons, total_ibs)
    
    def find_closest_ancients(
        self, 
        user_genotypes: Dict[str, str], 
        top_n: int = 20,
        min_snps: int = 5
    ) -> List[AncientMatch]:
        """
        Find the N closest ancient individuals to the user.
        
        Args:
            user_genotypes: Dict of rsid -> genotype
            top_n: Number of top matches to return
            min_snps: Minimum shared SNPs for a valid match
            
        Returns:
            List of AncientMatch objects, sorted by similarity (highest first)
        """
        matches = []
        
        for individual in self.individuals:
            ancient_snps = individual.get("snps", {})
            
            distance, shared, _ = self.calculate_ibs_distance(
                user_genotypes, ancient_snps, min_snps
            )
            
            if distance is None:
                continue
            
            similarity = 100 * (1 - distance)
            
            # Find shared traits
            shared_traits = self._find_shared_traits(user_genotypes, ancient_snps)
            
            match = AncientMatch(
                id=individual["id"],
                name=individual["name"],
                similarity=round(similarity, 2),
                distance=round(distance, 4),
                shared_snps=shared,
                total_snps=len(ancient_snps),
                period=individual.get("period", "Unknown"),
                culture=individual.get("culture", "Unknown"),
                site=individual.get("site", "Unknown"),
                country=individual.get("country", "Unknown"),
                age_text=individual.get("age_text", "Unknown"),
                description=individual.get("description", ""),
                lat=individual.get("lat"),
                lon=individual.get("lon"),
                mt_haplogroup=individual.get("mt_haplogroup"),
                y_haplogroup=individual.get("y_haplogroup"),
                shared_traits=shared_traits,
                paper=individual.get("paper"),
                pmid=individual.get("pmid")
            )
            matches.append(match)
        
        # Sort by similarity (highest first)
        matches.sort(key=lambda m: m.similarity, reverse=True)
        
        return matches[:top_n]
    
    def _find_shared_traits(
        self, 
        user_genotypes: Dict[str, str], 
        ancient_genotypes: Dict[str, str]
    ) -> List[str]:
        """Find phenotypic traits likely shared between user and ancient."""
        shared_traits = []
        
        # Key trait SNPs
        trait_snps = {
            "rs12913832": {
                "GG": "Blue eyes",
                "GA": "Blue/green eyes",
                "AA": "Brown eyes"
            },
            "rs1426654": {
                "AA": "Light skin (SLC24A5)",
                "GA": "Intermediate skin",
                "GG": "Dark skin (ancestral)"
            },
            "rs16891982": {
                "GG": "Light skin (SLC45A2)",
                "GC": "Intermediate skin",
                "CC": "Darker skin"
            },
            "rs4988235": {
                "AA": "Lactase persistent",
                "GA": "Lactase persistent",
                "GG": "Lactose intolerant"
            },
            "rs1805007": {
                "TT": "Red hair risk",
                "CT": "Red hair carrier",
                "CC": "Non-carrier"
            }
        }
        
        for rs, traits in trait_snps.items():
            user_geno = user_genotypes.get(rs, "").upper()
            ancient_geno = ancient_genotypes.get(rs, "").upper()
            
            if not user_geno or not ancient_geno:
                continue
            
            # Normalize heterozygotes
            if len(user_geno) == 2:
                user_geno = "".join(sorted(user_geno))
            if len(ancient_geno) == 2:
                ancient_geno = "".join(sorted(ancient_geno))
            
            user_trait = traits.get(user_geno)
            ancient_trait = traits.get(ancient_geno)
            
            if user_trait and user_trait == ancient_trait:
                shared_traits.append(user_trait)
        
        return shared_traits
    
    def match_to_cultures(
        self, 
        user_genotypes: Dict[str, str],
        min_snps: int = 3
    ) -> List[CultureAffinity]:
        """
        Calculate affinity scores to ancient cultures.
        
        Uses average similarity to individuals from each culture,
        plus diagnostic SNP matching.
        
        Args:
            user_genotypes: Dict of rsid -> genotype
            min_snps: Minimum shared SNPs per individual
            
        Returns:
            List of CultureAffinity objects, sorted by affinity (highest first)
        """
        # Group individuals by culture
        culture_individuals = defaultdict(list)
        for ind in self.individuals:
            culture_id = ind.get("culture", "Unknown").lower().replace(" ", "_").replace("(", "").replace(")", "")
            culture_individuals[culture_id].append(ind)
        
        affinities = []
        
        for culture_id, culture_data in self.cultures.items():
            # Find individuals matching this culture (fuzzy match)
            matching_individuals = []
            culture_name_lower = culture_data["name"].lower()
            culture_abbrev = culture_data.get("abbreviation", "").lower()
            
            for ind in self.individuals:
                ind_culture = ind.get("culture", "").lower()
                if (culture_name_lower in ind_culture or 
                    culture_abbrev in ind_culture or
                    culture_id in ind_culture.replace(" ", "_").replace("(", "").replace(")", "")):
                    matching_individuals.append(ind)
            
            if not matching_individuals:
                continue
            
            # Calculate average similarity to individuals
            similarities = []
            shared_snp_counts = []
            
            for ind in matching_individuals:
                ancient_snps = ind.get("snps", {})
                distance, shared, _ = self.calculate_ibs_distance(
                    user_genotypes, ancient_snps, min_snps
                )
                if distance is not None:
                    similarities.append(100 * (1 - distance))
                    shared_snp_counts.append(shared)
            
            if not similarities:
                continue
            
            avg_similarity = sum(similarities) / len(similarities)
            avg_shared = sum(shared_snp_counts) / len(shared_snp_counts)
            
            # Check diagnostic SNPs
            diagnostic_matches = {}
            diag_snps = culture_data.get("diagnostic_snps", {})
            for rs, info in diag_snps.items():
                user_geno = user_genotypes.get(rs, "").upper()
                if user_geno:
                    ancestral = info.get("ancestral", "")
                    if ancestral in user_geno:
                        diagnostic_matches[rs] = f"Carries {ancestral} ({info.get('note', '')})"
            
            affinity = CultureAffinity(
                id=culture_id,
                name=culture_data["name"],
                affinity=round(avg_similarity, 2),
                period=culture_data.get("period", "Unknown"),
                regions=culture_data.get("regions", []),
                description=culture_data.get("description", ""),
                sample_count=len(similarities),
                avg_shared_snps=round(avg_shared, 1),
                diagnostic_matches=diagnostic_matches,
                map_color=culture_data.get("map_color", "#888888")
            )
            affinities.append(affinity)
        
        # Sort by affinity (highest first)
        affinities.sort(key=lambda a: a.affinity, reverse=True)
        
        return affinities
    
    def generate_detailed_comparison(
        self, 
        user_genotypes: Dict[str, str], 
        ancient_id: str
    ) -> Dict[str, Any]:
        """
        Generate detailed SNP-by-SNP comparison with specific ancient individual.
        
        Args:
            user_genotypes: Dict of rsid -> genotype
            ancient_id: ID of ancient individual
            
        Returns:
            Dict with detailed comparison data
        """
        # Find the individual
        individual = None
        for ind in self.individuals:
            if ind["id"] == ancient_id:
                individual = ind
                break
        
        if not individual:
            return {"error": f"Individual {ancient_id} not found"}
        
        ancient_snps = individual.get("snps", {})
        shared_snps = set(user_genotypes.keys()) & set(ancient_snps.keys())
        
        comparisons = []
        for rs in sorted(shared_snps):
            user_geno = user_genotypes[rs].upper()
            ancient_geno = ancient_snps[rs].upper()
            
            # Calculate match
            user_alleles = set(user_geno) if len(user_geno) >= 2 else set()
            ancient_alleles = set(ancient_geno) if len(ancient_geno) >= 2 else set()
            
            if user_geno == ancient_geno or (user_alleles == ancient_alleles):
                match_type = "identical"
            elif len(user_alleles & ancient_alleles) > 0:
                match_type = "partial"
            else:
                match_type = "different"
            
            comparisons.append({
                "rsid": rs,
                "user": user_geno,
                "ancient": ancient_geno,
                "match_type": match_type
            })
        
        distance, shared, total_ibs = self.calculate_ibs_distance(
            user_genotypes, ancient_snps
        )
        
        return {
            "individual": {
                "id": individual["id"],
                "name": individual["name"],
                "period": individual.get("period"),
                "culture": individual.get("culture"),
                "site": individual.get("site"),
                "age_text": individual.get("age_text"),
                "description": individual.get("description"),
                "paper": individual.get("paper"),
                "pmid": individual.get("pmid")
            },
            "summary": {
                "similarity": round(100 * (1 - distance), 2) if distance else None,
                "distance": round(distance, 4) if distance else None,
                "shared_snps": shared,
                "total_ibs_score": total_ibs,
                "max_possible_ibs": 2 * shared
            },
            "comparisons": comparisons,
            "identical_count": sum(1 for c in comparisons if c["match_type"] == "identical"),
            "partial_count": sum(1 for c in comparisons if c["match_type"] == "partial"),
            "different_count": sum(1 for c in comparisons if c["match_type"] == "different")
        }
    
    def generate_timeline_data(
        self, 
        matches: List[AncientMatch]
    ) -> List[Dict]:
        """
        Generate data for timeline visualization.
        
        Args:
            matches: List of AncientMatch objects
            
        Returns:
            List of timeline entries with time periods and match info
        """
        timeline = []
        
        for match in matches:
            # Parse age from age_text
            age_text = match.age_text
            
            # Extract approximate year
            year = None
            if "BCE" in age_text:
                # Extract number before BCE
                import re
                nums = re.findall(r'[\d,]+', age_text)
                if nums:
                    year = -int(nums[0].replace(',', ''))
            elif "CE" in age_text:
                nums = re.findall(r'[\d,]+', age_text)
                if nums:
                    year = int(nums[0].replace(',', ''))
            
            if year is not None:
                timeline.append({
                    "name": match.name,
                    "year": year,
                    "year_display": age_text,
                    "period": match.period,
                    "culture": match.culture,
                    "similarity": match.similarity,
                    "lat": match.lat,
                    "lon": match.lon
                })
        
        # Sort by year (oldest first)
        timeline.sort(key=lambda x: x["year"])
        
        return timeline
    
    def get_all_periods(self) -> List[Dict]:
        """Get list of all time periods for reference."""
        periods = [
            {"name": "Mesolithic", "start": -10000, "end": -5000, "color": "#2E86AB"},
            {"name": "Neolithic", "start": -5000, "end": -2500, "color": "#F18F01"},
            {"name": "Chalcolithic", "start": -3500, "end": -2500, "color": "#D4A373"},
            {"name": "Bronze Age", "start": -2500, "end": -800, "color": "#CD7F32"},
            {"name": "Iron Age", "start": -800, "end": 43, "color": "#4A7C59"},
            {"name": "Roman", "start": 43, "end": 410, "color": "#DC143C"},
            {"name": "Early Medieval", "start": 410, "end": 1066, "color": "#6247AA"}
        ]
        return periods


def analyze_ancient_matches(user_genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Main function to run full ancient DNA analysis.
    
    Args:
        user_genotypes: Dict of rsid -> genotype
        
    Returns:
        Dict with all analysis results
    """
    matcher = AncientDNAMatcher()
    
    # Get matches
    matches = matcher.find_closest_ancients(user_genotypes, top_n=30)
    
    # Get culture affinities
    affinities = matcher.match_to_cultures(user_genotypes)
    
    # Generate timeline
    timeline = matcher.generate_timeline_data(matches)
    
    # Get periods reference
    periods = matcher.get_all_periods()
    
    return {
        "matches": [m.to_dict() for m in matches],
        "culture_affinities": [
            {
                "id": a.id,
                "name": a.name,
                "affinity": a.affinity,
                "period": a.period,
                "regions": a.regions,
                "description": a.description,
                "sample_count": a.sample_count,
                "map_color": a.map_color
            }
            for a in affinities
        ],
        "timeline": timeline,
        "periods": periods,
        "methodology": {
            "method": "Identity-by-State (IBS)",
            "description": "We compare your DNA directly to published ancient genomes by counting shared alleles at each SNP position.",
            "similarity_meaning": "100% = identical at all compared SNPs; 50% = typical for distantly related individuals",
            "note": "Unlike proprietary services, our methodology is fully transparent and every ancient sample is cited with its publication.",
            "limitations": [
                "Similarity depends on which SNPs are available in both samples",
                "Ancient DNA often has missing data due to degradation",
                "High similarity does not mean direct descent",
                "Population-level averages may not represent individuals"
            ]
        }
    }


# Export functions for use by other modules
__all__ = [
    'AncientDNAMatcher',
    'AncientMatch', 
    'CultureAffinity',
    'analyze_ancient_matches'
]


if __name__ == "__main__":
    # Test with sample genotypes
    sample_user = {
        "rs12913832": "GA",  # Blue/green eyes
        "rs1426654": "AA",   # Light skin
        "rs16891982": "GG",  # Light skin  
        "rs4988235": "GA",   # Lactase persistent
        "rs1805007": "CC",   # Non-red hair
        "rs1800407": "CT",
        "rs4778241": "GA",
        "rs12203592": "CC",
        "rs1042602": "CA"
    }
    
    results = analyze_ancient_matches(sample_user)
    
    print("=" * 70)
    print("ANCIENT DNA MATCHING RESULTS")
    print("=" * 70)
    
    print("\nTOP 10 CLOSEST ANCIENT INDIVIDUALS:")
    print("-" * 70)
    for i, match in enumerate(results["matches"][:10], 1):
        print(f"{i}. {match['name']} ({match['period']})")
        print(f"   Similarity: {match['similarity']}% | Shared SNPs: {match['shared_snps']}")
        print(f"   Culture: {match['culture']}")
        print(f"   Site: {match['site']}, {match['country']}")
        if match['shared_traits']:
            print(f"   Shared traits: {', '.join(match['shared_traits'])}")
        print()
    
    print("\nCULTURE AFFINITIES:")
    print("-" * 70)
    for aff in results["culture_affinities"][:10]:
        print(f"{aff['name']}: {aff['affinity']}% ({aff['period']})")
