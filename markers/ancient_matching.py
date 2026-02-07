"""
Ancient Individual Matching Module

Free "YourTrueAncestry" clone - matches your DNA against published ancient genomes.

This module provides:
- Genetic distance calculation using IBS (Identity by State)
- Matching against 30+ ancient individuals from published studies
- Cultural affinity scoring
- Detailed reports with historical context

Key principle: TRANSPARENCY. We show exactly how we calculate matches,
provide PMIDs for every ancient sample, and explain limitations clearly.

Sources:
- Reich Lab publications
- Allen Ancient DNA Resource (AADR)
- Published ancient DNA studies with PMID citations
"""

from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
from collections import defaultdict
import json
import math

# =============================================================================
# LOAD REFERENCE DATA
# =============================================================================

def load_ancient_individuals() -> Dict[str, Any]:
    """Load the ancient individuals database."""
    ref_path = Path(__file__).parent.parent / "references" / "ancient_individuals.json"
    
    if not ref_path.exists():
        return {"_metadata": {}}
    
    with open(ref_path, 'r') as f:
        return json.load(f)


def load_ancient_cultures() -> Dict[str, Any]:
    """Load the ancient cultures database."""
    ref_path = Path(__file__).parent.parent / "references" / "ancient_cultures.json"
    
    if not ref_path.exists():
        return {"_metadata": {}}
    
    with open(ref_path, 'r') as f:
        return json.load(f)


# Cache the data
_ANCIENT_INDIVIDUALS = None
_ANCIENT_CULTURES = None


def get_ancient_individuals() -> Dict[str, Any]:
    """Get cached ancient individuals data."""
    global _ANCIENT_INDIVIDUALS
    if _ANCIENT_INDIVIDUALS is None:
        _ANCIENT_INDIVIDUALS = load_ancient_individuals()
    return _ANCIENT_INDIVIDUALS


def get_ancient_cultures() -> Dict[str, Any]:
    """Get cached ancient cultures data."""
    global _ANCIENT_CULTURES
    if _ANCIENT_CULTURES is None:
        _ANCIENT_CULTURES = load_ancient_cultures()
    return _ANCIENT_CULTURES


# =============================================================================
# GENETIC DISTANCE CALCULATION
# =============================================================================

def normalize_genotype(genotype: str) -> str:
    """Normalize genotype to alphabetical order for comparison."""
    if not genotype:
        return ""
    if len(genotype) == 2:
        return "".join(sorted(genotype.upper()))
    return genotype.upper()


def calculate_genetic_distance(
    user_genos: Dict[str, str],
    ancient_genos: Dict[str, str],
    weighted: bool = True
) -> Optional[Dict[str, Any]]:
    """
    Calculate genetic distance between user and ancient individual using IBS.
    
    Identity by State (IBS) is a simple measure of genetic similarity:
    - 2 = identical genotype (e.g., AA vs AA)
    - 1 = partial match (e.g., AA vs AG)  
    - 0 = no match (e.g., AA vs GG)
    
    Args:
        user_genos: Dict mapping rsid -> genotype for user
        ancient_genos: Dict mapping rsid -> genotype for ancient individual
        weighted: Whether to weight informative SNPs more heavily
        
    Returns:
        Dict with distance, shared_snps, ibs_score, and details
        Returns None if insufficient shared SNPs
    """
    # Find shared SNPs (excluding null/missing in ancient)
    shared_snps = []
    for rsid in user_genos:
        if rsid in ancient_genos and ancient_genos[rsid]:
            user_geno = normalize_genotype(user_genos[rsid])
            ancient_geno = normalize_genotype(ancient_genos[rsid])
            if user_geno and ancient_geno:
                shared_snps.append(rsid)
    
    if len(shared_snps) < 3:
        return None
    
    # Calculate IBS
    total_ibs = 0
    max_ibs = 0
    details = []
    
    # SNP weights (ancestry-informative markers weighted higher)
    informative_snps = {
        "rs1426654": 2.0,   # SLC24A5 - highly informative
        "rs16891982": 2.0,  # SLC45A2 - highly informative
        "rs12913832": 1.5,  # HERC2 - eye color
        "rs4988235": 1.5,   # LCT - lactase
        "rs3827760": 2.0,   # EDAR - East Asian specific
        "rs2814778": 2.0,   # DARC - African specific
    }
    
    for rsid in shared_snps:
        user_geno = normalize_genotype(user_genos[rsid])
        ancient_geno = normalize_genotype(ancient_genos[rsid])
        
        weight = informative_snps.get(rsid, 1.0) if weighted else 1.0
        
        # Calculate IBS (0, 1, or 2)
        if user_geno == ancient_geno:
            ibs = 2
        elif set(user_geno) & set(ancient_geno):  # At least one allele shared
            ibs = 1
        else:
            ibs = 0
        
        total_ibs += ibs * weight
        max_ibs += 2 * weight
        
        details.append({
            "rsid": rsid,
            "user": user_geno,
            "ancient": ancient_geno,
            "ibs": ibs,
            "match": "identical" if ibs == 2 else "partial" if ibs == 1 else "different"
        })
    
    # Calculate normalized distance (0 = identical, 1 = completely different)
    ibs_score = total_ibs / max_ibs if max_ibs > 0 else 0
    distance = 1 - ibs_score
    
    return {
        "distance": round(distance, 4),
        "similarity": round(ibs_score * 100, 1),  # As percentage
        "shared_snps": len(shared_snps),
        "total_ibs": total_ibs,
        "max_ibs": max_ibs,
        "details": details
    }


def calculate_trait_matches(
    user_genos: Dict[str, str],
    ancient_data: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """
    Determine which physical traits user shares with ancient individual.
    
    Returns list of matched/different traits based on known marker effects.
    """
    trait_markers = {
        "rs1426654": {
            "trait": "Skin pigmentation",
            "derived": "AA",
            "derived_effect": "Light skin",
            "ancestral": "GG",
            "ancestral_effect": "Dark skin"
        },
        "rs16891982": {
            "trait": "Skin pigmentation",
            "derived": "GG",
            "derived_effect": "Light skin",
            "ancestral": "CC",
            "ancestral_effect": "Darker skin"
        },
        "rs12913832": {
            "trait": "Eye color",
            "derived": "GG",
            "derived_effect": "Blue eyes",
            "ancestral": "AA",
            "ancestral_effect": "Brown eyes"
        },
        "rs4988235": {
            "trait": "Lactase persistence",
            "derived": "AA",
            "derived_effect": "Can digest milk as adult",
            "ancestral": "GG",
            "ancestral_effect": "Lactose intolerant"
        }
    }
    
    ancient_snps = ancient_data.get("snps", {})
    matches = []
    
    for rsid, info in trait_markers.items():
        user_geno = normalize_genotype(user_genos.get(rsid, ""))
        ancient_geno = normalize_genotype(ancient_snps.get(rsid, ""))
        
        if not user_geno or not ancient_geno:
            continue
        
        # Determine user's trait
        if user_geno == info["derived"]:
            user_trait = info["derived_effect"]
        elif user_geno == info["ancestral"]:
            user_trait = info["ancestral_effect"]
        else:
            user_trait = "Intermediate"
        
        # Determine ancient's trait
        if ancient_geno == info["derived"]:
            ancient_trait = info["derived_effect"]
        elif ancient_geno == info["ancestral"]:
            ancient_trait = info["ancestral_effect"]
        else:
            ancient_trait = "Intermediate"
        
        match = user_geno == ancient_geno
        
        matches.append({
            "trait": info["trait"],
            "rsid": rsid,
            "user_genotype": user_geno,
            "user_trait": user_trait,
            "ancient_genotype": ancient_geno,
            "ancient_trait": ancient_trait,
            "match": match,
            "symbol": "✓" if match else "✗"
        })
    
    return matches


# =============================================================================
# ANCIENT INDIVIDUAL MATCHING
# =============================================================================

def find_closest_ancients(
    user_genos: Dict[str, str],
    top_n: int = 10,
    min_shared_snps: int = 3
) -> List[Dict[str, Any]]:
    """
    Find the closest ancient individuals to the user.
    
    Args:
        user_genos: User's genotypes
        top_n: Number of top matches to return
        min_shared_snps: Minimum SNPs required for a valid comparison
        
    Returns:
        List of ancient matches sorted by similarity
    """
    ancient_data = get_ancient_individuals()
    results = []
    
    for ancient_id, individual in ancient_data.items():
        if ancient_id.startswith("_"):  # Skip metadata
            continue
        
        ancient_snps = individual.get("snps", {})
        if not ancient_snps:
            continue
        
        distance_result = calculate_genetic_distance(user_genos, ancient_snps)
        
        if distance_result is None:
            continue
        
        if distance_result["shared_snps"] < min_shared_snps:
            continue
        
        trait_matches = calculate_trait_matches(user_genos, individual)
        
        results.append({
            "id": ancient_id,
            "name": individual.get("name", ancient_id),
            "site": individual.get("site", "Unknown"),
            "country": individual.get("country", ""),
            "age_bp": individual.get("age_bp", 0),
            "age_calBCE": individual.get("age_calBCE", 0),
            "period": individual.get("period", "Unknown"),
            "culture": individual.get("culture", ""),
            "culture_name": individual.get("culture_name", ""),
            "sex": individual.get("sex", "Unknown"),
            "description": individual.get("description", ""),
            "significance": individual.get("significance", ""),
            "physical_traits": individual.get("physical_traits", []),
            "paper": individual.get("paper", ""),
            "pmid": individual.get("pmid", ""),
            
            # Matching results
            "genetic_distance": distance_result["distance"],
            "similarity_percent": distance_result["similarity"],
            "shared_snps": distance_result["shared_snps"],
            "trait_comparison": trait_matches,
            "snp_details": distance_result["details"]
        })
    
    # Sort by similarity (highest first)
    results.sort(key=lambda x: x["similarity_percent"], reverse=True)
    
    return results[:top_n]


# =============================================================================
# CULTURAL AFFINITY
# =============================================================================

def match_to_cultures(user_genos: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Calculate affinity scores for each ancient culture.
    
    Aggregates individual matches within each culture to give overall affinity.
    
    Args:
        user_genos: User's genotypes
        
    Returns:
        List of cultures with affinity scores
    """
    ancient_individuals = get_ancient_individuals()
    ancient_cultures = get_ancient_cultures()
    
    # Group individuals by culture and calculate average similarity
    culture_scores = defaultdict(lambda: {"total_similarity": 0, "count": 0, "individuals": []})
    
    for ancient_id, individual in ancient_individuals.items():
        if ancient_id.startswith("_"):
            continue
        
        culture = individual.get("culture", "")
        if not culture:
            continue
        
        ancient_snps = individual.get("snps", {})
        distance_result = calculate_genetic_distance(user_genos, ancient_snps)
        
        if distance_result is None:
            continue
        
        culture_scores[culture]["total_similarity"] += distance_result["similarity"]
        culture_scores[culture]["count"] += 1
        culture_scores[culture]["individuals"].append({
            "id": ancient_id,
            "name": individual.get("name", ancient_id),
            "similarity": distance_result["similarity"]
        })
    
    # Build results with culture metadata
    results = []
    
    for culture_id, scores in culture_scores.items():
        if scores["count"] == 0:
            continue
        
        avg_similarity = scores["total_similarity"] / scores["count"]
        culture_info = ancient_cultures.get(culture_id, {})
        
        results.append({
            "culture_id": culture_id,
            "name": culture_info.get("name", culture_id),
            "full_name": culture_info.get("full_name", culture_id),
            "period_start": culture_info.get("period_start", 0),
            "period_end": culture_info.get("period_end", 0),
            "color": culture_info.get("color", "#888888"),
            "icon": culture_info.get("icon", "🏛️"),
            "regions": culture_info.get("regions", []),
            "description": culture_info.get("description", ""),
            "modern_contribution": culture_info.get("modern_contribution", ""),
            "pmid": culture_info.get("pmid", []),
            
            # Affinity results
            "affinity_percent": round(avg_similarity, 1),
            "individuals_compared": scores["count"],
            "top_individuals": sorted(scores["individuals"], key=lambda x: x["similarity"], reverse=True)[:3]
        })
    
    # Sort by affinity
    results.sort(key=lambda x: x["affinity_percent"], reverse=True)
    
    return results


# =============================================================================
# REPORT GENERATION
# =============================================================================

def make_bar(pct: float, width: int = 20) -> str:
    """Create a simple text bar chart."""
    filled = int(pct / 100 * width)
    return "█" * filled + "░" * (width - filled)


def generate_ancient_matches_report(user_genos: Dict[str, str]) -> str:
    """
    Generate a detailed text report of ancient matches.
    
    Args:
        user_genos: User's genotypes
        
    Returns:
        Formatted text report
    """
    lines = []
    lines.append("=" * 70)
    lines.append("🏛️ YOUR ANCIENT DNA MATCHES")
    lines.append("=" * 70)
    lines.append("")
    lines.append("This analysis compares your DNA against published ancient genomes")
    lines.append("using Identity by State (IBS) - a measure of shared alleles.")
    lines.append("")
    lines.append("⚠️  IMPORTANT: This is an EDUCATIONAL tool, not precise ancestry.")
    lines.append("Ancient DNA coverage is sparse - we only compare overlapping SNPs.")
    lines.append("")
    
    # Get matches
    matches = find_closest_ancients(user_genos, top_n=10)
    
    if not matches:
        lines.append("❌ Insufficient SNP overlap with ancient reference data.")
        return "\n".join(lines)
    
    lines.append("-" * 70)
    lines.append("🏆 TOP ANCIENT MATCHES")
    lines.append("-" * 70)
    lines.append("")
    
    for i, match in enumerate(matches[:5], 1):
        lines.append(f"#{i} {match['name'].upper()}")
        lines.append(f"    Genetic Distance: {match['genetic_distance']:.2f}")
        lines.append(f"    Similarity: {match['similarity_percent']:.1f}%  {make_bar(match['similarity_percent'])}")
        lines.append(f"    Period: {match['period']} ({match['age_calBCE']} BCE)")
        lines.append(f"    Site: {match['site']}")
        lines.append(f"    Culture: {match['culture_name']}")
        
        if match['description']:
            # Word wrap description
            desc = match['description']
            words = desc.split()
            line = "    "
            for word in words:
                if len(line) + len(word) > 65:
                    lines.append(line)
                    line = "    "
                line += word + " "
            if line.strip():
                lines.append(line)
        
        # Trait comparison
        traits = [t for t in match.get('trait_comparison', []) if t['match']]
        if traits:
            lines.append("\n    Shared traits:")
            for t in traits:
                lines.append(f"      ✓ {t['trait']}: {t['user_trait']}")
        
        if match.get('pmid'):
            lines.append(f"\n    Reference: PMID:{match['pmid']}")
        
        lines.append("")
    
    # Cultural affinities
    lines.append("-" * 70)
    lines.append("📊 CULTURAL AFFINITIES")
    lines.append("-" * 70)
    lines.append("")
    
    cultures = match_to_cultures(user_genos)
    
    for culture in cultures[:8]:
        name = culture['name']
        pct = culture['affinity_percent']
        icon = culture.get('icon', '🏛️')
        bar = make_bar(pct)
        lines.append(f"{icon} {name:30} {pct:5.1f}%  {bar}")
    
    lines.append("")
    
    # Timeline
    lines.append("-" * 70)
    lines.append("📅 YOUR MATCHES THROUGH TIME")
    lines.append("-" * 70)
    lines.append("")
    
    # Sort matches by age
    sorted_matches = sorted(matches, key=lambda x: x.get('age_calBCE', 0), reverse=True)
    
    for match in sorted_matches[:7]:
        year = match.get('age_calBCE', 0)
        year_str = f"{abs(year)} BCE" if year < 0 else f"{year} CE"
        if year == 0:
            year_str = "Unknown"
        lines.append(f"  {year_str:12} │ {match['name']} ({match['culture_name']})")
    
    lines.append("")
    
    # Methodology
    lines.append("=" * 70)
    lines.append("📖 METHODOLOGY")
    lines.append("=" * 70)
    lines.append("")
    lines.append("How we calculate genetic distance:")
    lines.append("  1. Find SNPs present in both your data and ancient sample")
    lines.append("  2. Compare genotypes using Identity by State (IBS):")
    lines.append("     - IBS=2: Identical genotype (AA vs AA)")
    lines.append("     - IBS=1: Partial match (AA vs AG)")
    lines.append("     - IBS=0: No match (AA vs GG)")
    lines.append("  3. Sum IBS scores and normalize to 0-1 scale")
    lines.append("  4. Distance = 1 - normalized_IBS")
    lines.append("")
    lines.append("Ancestry-informative SNPs (SLC24A5, SLC45A2, HERC2, LCT)")
    lines.append("are weighted more heavily as they better distinguish populations.")
    lines.append("")
    lines.append("LIMITATIONS:")
    lines.append("• Consumer arrays only test ~1M of 3B SNPs")
    lines.append("• Ancient DNA often has missing data")
    lines.append("• We only compare 5-20 overlapping SNPs per sample")
    lines.append("• This is educational, not definitive ancestry")
    lines.append("")
    lines.append("All ancient samples from peer-reviewed publications.")
    lines.append("PMIDs provided for verification.")
    lines.append("")
    
    return "\n".join(lines)


# =============================================================================
# JSON EXPORT FOR DASHBOARD
# =============================================================================

def get_ancient_matches_json(user_genos: Dict[str, str]) -> Dict[str, Any]:
    """
    Get comprehensive ancient matching data as JSON for dashboard.
    
    Args:
        user_genos: User's genotypes
        
    Returns:
        Structured data for dashboard visualization
    """
    matches = find_closest_ancients(user_genos, top_n=15)
    cultures = match_to_cultures(user_genos)
    
    # Build timeline data
    timeline = []
    for match in matches:
        year = match.get('age_calBCE', 0)
        timeline.append({
            "year": year,
            "year_display": f"{abs(year)} BCE" if year < 0 else f"{year} CE",
            "name": match['name'],
            "culture": match['culture_name'],
            "similarity": match['similarity_percent']
        })
    timeline.sort(key=lambda x: x['year'], reverse=True)
    
    # Get culture metadata for visualization
    culture_data = get_ancient_cultures()
    
    return {
        "top_matches": matches[:10],
        "all_matches": matches,
        "culture_affinities": cultures,
        "timeline": timeline,
        "culture_metadata": {k: v for k, v in culture_data.items() if not k.startswith("_")},
        "statistics": {
            "total_individuals_compared": len(matches),
            "total_cultures_compared": len(cultures),
            "best_match_similarity": matches[0]["similarity_percent"] if matches else 0,
            "avg_snps_compared": sum(m["shared_snps"] for m in matches) / len(matches) if matches else 0
        },
        "methodology": {
            "method": "Identity by State (IBS)",
            "description": "Counts shared alleles between user and ancient samples",
            "limitations": [
                "Consumer arrays test ~0.03% of the genome",
                "Ancient DNA has missing data",
                "Only 5-20 SNPs typically overlap",
                "Educational tool, not precise ancestry"
            ],
            "sources": [
                "Allen Ancient DNA Resource (AADR)",
                "Reich Lab publications",
                "Published studies with PMID citations"
            ]
        }
    }


# =============================================================================
# CONVENIENCE EXPORTS
# =============================================================================

def analyze_ancient_ancestry(user_genos: Dict[str, str]) -> Dict[str, Any]:
    """
    Complete ancient ancestry analysis for integration with main analysis.
    
    Args:
        user_genos: User's genotypes
        
    Returns:
        Combined results suitable for agent_summary.json
    """
    json_data = get_ancient_matches_json(user_genos)
    
    # Simplify for agent summary
    top_3 = []
    for m in json_data.get("top_matches", [])[:3]:
        top_3.append({
            "name": m["name"],
            "culture": m["culture_name"],
            "period": m["period"],
            "similarity": m["similarity_percent"],
            "pmid": m.get("pmid", "")
        })
    
    top_cultures = []
    for c in json_data.get("culture_affinities", [])[:5]:
        top_cultures.append({
            "name": c["name"],
            "affinity": c["affinity_percent"],
            "icon": c.get("icon", "🏛️")
        })
    
    return {
        "top_ancient_matches": top_3,
        "cultural_affinities": top_cultures,
        "total_compared": json_data.get("statistics", {}).get("total_individuals_compared", 0),
        "methodology_note": "IBS-based genetic distance. Educational, not precise ancestry.",
        "full_data": json_data
    }
