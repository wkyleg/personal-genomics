"""
Ancient Individual Matching Module - v2.0

Free "YourTrueAncestry" clone - matches your DNA against published ancient genomes.

MAJOR UPDATE v2.0:
- Expanded SNP coverage (30-40 SNPs per ancient individual)
- Confidence levels based on shared SNP count
- Meaningful similarity display with N (sample size)
- Percentile rankings across all ancients

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
# CONFIGURATION
# =============================================================================

# Confidence thresholds for shared SNPs
CONFIDENCE_HIGH = 30      # High confidence: 30+ shared SNPs
CONFIDENCE_MEDIUM = 20    # Medium confidence: 20-29 shared SNPs
CONFIDENCE_LOW = 10       # Low confidence: 10-19 shared SNPs
MINIMUM_SNPS = 5          # Below this, don't report match

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
    # Handle special cases like "delC" for deletion
    if "del" in genotype.lower():
        return genotype
    if len(genotype) == 2:
        return "".join(sorted(genotype.upper()))
    return genotype.upper()


def get_confidence_level(shared_snps: int) -> str:
    """
    Determine confidence level based on number of shared SNPs.
    
    Returns: 'high', 'medium', 'low', or 'insufficient'
    """
    if shared_snps >= CONFIDENCE_HIGH:
        return "high"
    elif shared_snps >= CONFIDENCE_MEDIUM:
        return "medium"
    elif shared_snps >= CONFIDENCE_LOW:
        return "low"
    elif shared_snps >= MINIMUM_SNPS:
        return "very_low"
    else:
        return "insufficient"


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
        Dict with distance, shared_snps, ibs_score, confidence, and details
        Returns None if insufficient shared SNPs
    """
    # Find shared SNPs (excluding null/missing in ancient)
    shared_snps = []
    for rsid in user_genos:
        if rsid in ancient_genos and ancient_genos[rsid]:
            user_geno = normalize_genotype(user_genos[rsid])
            ancient_geno = normalize_genotype(ancient_genos[rsid])
            if user_geno and ancient_geno and "del" not in ancient_geno.lower():
                shared_snps.append(rsid)
    
    if len(shared_snps) < MINIMUM_SNPS:
        return None
    
    # Calculate IBS
    total_ibs = 0
    max_ibs = 0
    identical_count = 0
    partial_count = 0
    different_count = 0
    details = []
    
    # SNP weights (ancestry-informative markers weighted higher)
    informative_snps = {
        # Pigmentation - highly ancestry-informative
        "rs1426654": 2.0,   # SLC24A5 - highly informative
        "rs16891982": 2.0,  # SLC45A2 - highly informative
        "rs12913832": 1.5,  # HERC2 - eye color
        "rs1042602": 1.3,   # TYR - pigmentation
        "rs1800407": 1.3,   # OCA2 - eyes
        "rs7495174": 1.3,   # OCA2 - eyes
        # Diet/metabolism
        "rs4988235": 1.5,   # LCT - lactase
        "rs174546": 1.3,    # FADS1 - fatty acids
        # Population-specific
        "rs3827760": 2.0,   # EDAR - East Asian specific
        "rs2814778": 2.0,   # DARC - African specific
        "rs1800414": 1.5,   # OCA2 - East Asian
    }
    
    for rsid in shared_snps:
        user_geno = normalize_genotype(user_genos[rsid])
        ancient_geno = normalize_genotype(ancient_genos[rsid])
        
        weight = informative_snps.get(rsid, 1.0) if weighted else 1.0
        
        # Calculate IBS (0, 1, or 2)
        if user_geno == ancient_geno:
            ibs = 2
            identical_count += 1
        elif set(user_geno) & set(ancient_geno):  # At least one allele shared
            ibs = 1
            partial_count += 1
        else:
            ibs = 0
            different_count += 1
        
        total_ibs += ibs * weight
        max_ibs += 2 * weight
        
        details.append({
            "rsid": rsid,
            "user": user_geno,
            "ancient": ancient_geno,
            "ibs": ibs,
            "weight": weight,
            "match": "identical" if ibs == 2 else "partial" if ibs == 1 else "different"
        })
    
    # Calculate normalized distance (0 = identical, 1 = completely different)
    ibs_score = total_ibs / max_ibs if max_ibs > 0 else 0
    distance = 1 - ibs_score
    
    # Determine confidence level
    confidence = get_confidence_level(len(shared_snps))
    
    return {
        "distance": round(distance, 4),
        "similarity": round(ibs_score * 100, 1),  # As percentage
        "shared_snps": len(shared_snps),
        "identical_snps": identical_count,
        "partial_snps": partial_count,
        "different_snps": different_count,
        "total_ibs": total_ibs,
        "max_ibs": max_ibs,
        "confidence": confidence,
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
            "trait": "Skin pigmentation (SLC24A5)",
            "derived": "AA",
            "derived_effect": "Light skin",
            "ancestral": "GG",
            "ancestral_effect": "Dark skin"
        },
        "rs16891982": {
            "trait": "Skin pigmentation (SLC45A2)",
            "derived": "GG",
            "derived_effect": "Light skin",
            "ancestral": "CC",
            "ancestral_effect": "Darker skin"
        },
        "rs12913832": {
            "trait": "Eye color (HERC2)",
            "derived": "GG",
            "derived_effect": "Blue eyes",
            "ancestral": "AA",
            "ancestral_effect": "Brown eyes"
        },
        "rs4988235": {
            "trait": "Lactase persistence (LCT)",
            "derived": "AA",
            "derived_effect": "Can digest milk as adult",
            "ancestral": "GG",
            "ancestral_effect": "Lactose intolerant"
        },
        "rs1805007": {
            "trait": "Hair/skin (MC1R R151C)",
            "derived": "TT",
            "derived_effect": "Red hair/fair skin risk",
            "ancestral": "CC",
            "ancestral_effect": "Normal pigmentation"
        },
        "rs1805008": {
            "trait": "Hair/skin (MC1R R160W)",
            "derived": "TT",
            "derived_effect": "Red hair/fair skin risk",
            "ancestral": "CC",
            "ancestral_effect": "Normal pigmentation"
        },
        "rs17822931": {
            "trait": "Earwax type (ABCC11)",
            "derived": "CC",
            "derived_effect": "Dry earwax (East Asian)",
            "ancestral": "TT",
            "ancestral_effect": "Wet earwax (European/African)"
        },
        "rs3827760": {
            "trait": "Hair thickness (EDAR)",
            "derived": "GG",
            "derived_effect": "Thick hair (East Asian)",
            "ancestral": "AA",
            "ancestral_effect": "Normal hair"
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
            user_trait = "Intermediate/heterozygous"
        
        # Determine ancient's trait
        if ancient_geno == info["derived"]:
            ancient_trait = info["derived_effect"]
        elif ancient_geno == info["ancestral"]:
            ancient_trait = info["ancestral_effect"]
        else:
            ancient_trait = "Intermediate/heterozygous"
        
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
    top_n: int = 15,
    min_confidence: str = "very_low"
) -> List[Dict[str, Any]]:
    """
    Find the closest ancient individuals to the user.
    
    Args:
        user_genos: User's genotypes
        top_n: Number of top matches to return
        min_confidence: Minimum confidence level ('high', 'medium', 'low', 'very_low')
        
    Returns:
        List of ancient matches sorted by similarity, with percentile rankings
    """
    ancient_data = get_ancient_individuals()
    results = []
    
    # Confidence level ordering
    confidence_order = {"insufficient": 0, "very_low": 1, "low": 2, "medium": 3, "high": 4}
    min_conf_value = confidence_order.get(min_confidence, 1)
    
    for ancient_id, individual in ancient_data.items():
        if ancient_id.startswith("_"):  # Skip metadata
            continue
        
        ancient_snps = individual.get("snps", {})
        if not ancient_snps:
            continue
        
        distance_result = calculate_genetic_distance(user_genos, ancient_snps)
        
        if distance_result is None:
            continue
        
        # Filter by confidence
        conf_value = confidence_order.get(distance_result["confidence"], 0)
        if conf_value < min_conf_value:
            continue
        
        trait_matches = calculate_trait_matches(user_genos, individual)
        
        results.append({
            "id": ancient_id,
            "name": individual.get("name", ancient_id),
            "site": individual.get("site", "Unknown"),
            "country": individual.get("country", ""),
            "coordinates": individual.get("coordinates", {}),
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
            "coverage": individual.get("coverage", "unknown"),
            
            # Matching results - NEW: more detailed metrics
            "genetic_distance": distance_result["distance"],
            "similarity_percent": distance_result["similarity"],
            "shared_snps": distance_result["shared_snps"],
            "identical_snps": distance_result["identical_snps"],
            "partial_snps": distance_result["partial_snps"],
            "different_snps": distance_result["different_snps"],
            "confidence": distance_result["confidence"],
            "trait_comparison": trait_matches,
            "snp_details": distance_result["details"]
        })
    
    # Sort by similarity (highest first)
    results.sort(key=lambda x: x["similarity_percent"], reverse=True)
    
    # Calculate percentile rankings
    total_matches = len(results)
    for i, result in enumerate(results):
        # Percentile: how many are BELOW this individual
        percentile = ((total_matches - i - 1) / total_matches) * 100 if total_matches > 1 else 50
        result["rank"] = i + 1
        result["total_compared"] = total_matches
        result["percentile"] = round(percentile, 1)
        result["rank_display"] = f"#{i + 1} of {total_matches}"
    
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
    culture_scores = defaultdict(lambda: {
        "total_similarity": 0, 
        "count": 0, 
        "individuals": [],
        "total_shared_snps": 0,
        "min_shared": float('inf'),
        "max_shared": 0
    })
    
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
        
        scores = culture_scores[culture]
        scores["total_similarity"] += distance_result["similarity"]
        scores["count"] += 1
        scores["total_shared_snps"] += distance_result["shared_snps"]
        scores["min_shared"] = min(scores["min_shared"], distance_result["shared_snps"])
        scores["max_shared"] = max(scores["max_shared"], distance_result["shared_snps"])
        scores["individuals"].append({
            "id": ancient_id,
            "name": individual.get("name", ancient_id),
            "similarity": distance_result["similarity"],
            "shared_snps": distance_result["shared_snps"],
            "confidence": distance_result["confidence"]
        })
    
    # Build results with culture metadata
    results = []
    
    for culture_id, scores in culture_scores.items():
        if scores["count"] == 0:
            continue
        
        avg_similarity = scores["total_similarity"] / scores["count"]
        avg_shared = scores["total_shared_snps"] / scores["count"]
        culture_info = ancient_cultures.get(culture_id, {})
        
        # Determine overall confidence for culture
        if scores["min_shared"] >= CONFIDENCE_HIGH:
            confidence = "high"
        elif scores["min_shared"] >= CONFIDENCE_MEDIUM:
            confidence = "medium"
        elif scores["min_shared"] >= CONFIDENCE_LOW:
            confidence = "low"
        else:
            confidence = "very_low"
        
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
            
            # Affinity results - NEW: more detailed
            "affinity_percent": round(avg_similarity, 1),
            "individuals_compared": scores["count"],
            "avg_shared_snps": round(avg_shared, 1),
            "shared_snps_range": f"{scores['min_shared']}-{scores['max_shared']}",
            "confidence": confidence,
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


def confidence_emoji(confidence: str) -> str:
    """Return emoji for confidence level."""
    return {
        "high": "🟢",
        "medium": "🟡", 
        "low": "🟠",
        "very_low": "🔴",
        "insufficient": "⚫"
    }.get(confidence, "⚪")


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
    lines.append("Results depend on the number of shared SNPs between your data and")
    lines.append("each ancient sample. We show N (sample size) for transparency.")
    lines.append("")
    lines.append("CONFIDENCE LEVELS:")
    lines.append("  🟢 High (30+ SNPs)    🟡 Medium (20-29)    🟠 Low (10-19)")
    lines.append("")
    
    # Get matches
    matches = find_closest_ancients(user_genos, top_n=15)
    
    if not matches:
        lines.append("❌ Insufficient SNP overlap with ancient reference data.")
        return "\n".join(lines)
    
    lines.append("-" * 70)
    lines.append("🏆 TOP ANCIENT MATCHES (with statistical details)")
    lines.append("-" * 70)
    lines.append("")
    
    for i, match in enumerate(matches[:7], 1):
        conf = confidence_emoji(match['confidence'])
        lines.append(f"#{i} {match['name'].upper()}")
        lines.append(f"    {conf} Confidence: {match['confidence'].upper()} ({match['shared_snps']} shared SNPs)")
        lines.append(f"    Shared genotypes: {match['identical_snps']}/{match['shared_snps']} identical ({match['similarity_percent']:.1f}%)")
        lines.append(f"    {make_bar(match['similarity_percent'])}")
        lines.append(f"    Rank: {match['rank_display']}")
        lines.append(f"    Period: {match['period']} ({match['age_calBCE']} BCE)")
        lines.append(f"    Site: {match['site']}, {match['country']}")
        lines.append(f"    Culture: {match['culture_name']}")
        
        if match['description']:
            # Word wrap description
            desc = match['description'][:150]
            if len(match['description']) > 150:
                desc += "..."
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
            for t in traits[:4]:
                lines.append(f"      ✓ {t['trait']}: {t['user_trait']}")
        
        if match.get('pmid'):
            lines.append(f"\n    Reference: PMID:{match['pmid']}")
        
        lines.append("")
    
    # Cultural affinities
    lines.append("-" * 70)
    lines.append("📊 CULTURAL AFFINITIES (by population group)")
    lines.append("-" * 70)
    lines.append("")
    
    cultures = match_to_cultures(user_genos)
    
    for culture in cultures[:8]:
        name = culture['name']
        pct = culture['affinity_percent']
        icon = culture.get('icon', '🏛️')
        conf = confidence_emoji(culture['confidence'])
        snp_range = culture['shared_snps_range']
        bar = make_bar(pct)
        lines.append(f"{icon} {name:25} {pct:5.1f}%  {bar}  {conf} ({snp_range} SNPs)")
    
    lines.append("")
    
    # Timeline
    lines.append("-" * 70)
    lines.append("📅 YOUR MATCHES THROUGH TIME")
    lines.append("-" * 70)
    lines.append("")
    
    # Sort matches by age
    sorted_matches = sorted(matches, key=lambda x: x.get('age_calBCE', 0), reverse=True)
    
    for match in sorted_matches[:8]:
        year = match.get('age_calBCE', 0)
        year_str = f"{abs(year)} BCE" if year > 0 else f"{abs(year)} CE"
        if year == 0:
            year_str = "Unknown"
        conf = confidence_emoji(match['confidence'])
        lines.append(f"  {year_str:12} │ {match['name']} ({match['culture_name']}) {conf}")
    
    lines.append("")
    
    # Statistics summary
    lines.append("-" * 70)
    lines.append("📈 STATISTICS SUMMARY")
    lines.append("-" * 70)
    lines.append("")
    
    total_compared = matches[0].get('total_compared', len(matches)) if matches else 0
    avg_snps = sum(m['shared_snps'] for m in matches) / len(matches) if matches else 0
    high_conf = len([m for m in matches if m['confidence'] == 'high'])
    med_conf = len([m for m in matches if m['confidence'] == 'medium'])
    
    lines.append(f"  Total ancient individuals compared: {total_compared}")
    lines.append(f"  Average SNPs shared per individual: {avg_snps:.1f}")
    lines.append(f"  High-confidence matches (30+ SNPs): {high_conf}")
    lines.append(f"  Medium-confidence matches (20-29): {med_conf}")
    lines.append(f"  Best match similarity: {matches[0]['similarity_percent']:.1f}% ({matches[0]['shared_snps']} SNPs)")
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
    lines.append("  4. Similarity% = normalized_IBS × 100")
    lines.append("")
    lines.append("Ancestry-informative SNPs (SLC24A5, SLC45A2, HERC2, LCT)")
    lines.append("are weighted more heavily as they better distinguish populations.")
    lines.append("")
    lines.append("CONFIDENCE LEVELS based on shared SNP count:")
    lines.append("  🟢 HIGH (30+): Strong statistical power, reliable results")
    lines.append("  🟡 MEDIUM (20-29): Good for general patterns")
    lines.append("  🟠 LOW (10-19): Suggestive but interpret with caution")
    lines.append("  🔴 VERY LOW (<10): Limited reliability")
    lines.append("")
    lines.append("LIMITATIONS:")
    lines.append("• Consumer arrays only test ~0.03% of genome")
    lines.append("• Ancient DNA has missing/low-quality data")
    lines.append("• High similarity ≠ direct descent")
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
    matches = find_closest_ancients(user_genos, top_n=20)
    cultures = match_to_cultures(user_genos)
    
    # Build timeline data
    timeline = []
    for match in matches:
        year = match.get('age_calBCE', 0)
        timeline.append({
            "year": year,
            "year_display": f"{abs(year)} BCE" if year > 0 else f"{year} CE",
            "name": match['name'],
            "culture": match['culture_name'],
            "similarity": match['similarity_percent'],
            "shared_snps": match['shared_snps'],
            "confidence": match['confidence']
        })
    timeline.sort(key=lambda x: x['year'], reverse=True)
    
    # Get culture metadata for visualization
    culture_data = get_ancient_cultures()
    
    # Calculate overall statistics
    if matches:
        avg_snps = sum(m['shared_snps'] for m in matches) / len(matches)
        max_snps = max(m['shared_snps'] for m in matches)
        min_snps = min(m['shared_snps'] for m in matches)
        high_conf_count = len([m for m in matches if m['confidence'] == 'high'])
        med_conf_count = len([m for m in matches if m['confidence'] == 'medium'])
    else:
        avg_snps = max_snps = min_snps = 0
        high_conf_count = med_conf_count = 0
    
    return {
        "top_matches": matches[:10],
        "all_matches": matches,
        "culture_affinities": cultures,
        "timeline": timeline,
        "culture_metadata": {k: v for k, v in culture_data.items() if not k.startswith("_")},
        "statistics": {
            "total_individuals_compared": matches[0].get('total_compared', len(matches)) if matches else 0,
            "total_cultures_compared": len(cultures),
            "best_match_similarity": matches[0]["similarity_percent"] if matches else 0,
            "best_match_name": matches[0]["name"] if matches else "N/A",
            "avg_snps_compared": round(avg_snps, 1),
            "snps_range": f"{min_snps}-{max_snps}",
            "high_confidence_matches": high_conf_count,
            "medium_confidence_matches": med_conf_count
        },
        "confidence_thresholds": {
            "high": CONFIDENCE_HIGH,
            "medium": CONFIDENCE_MEDIUM,
            "low": CONFIDENCE_LOW,
            "minimum": MINIMUM_SNPS
        },
        "methodology": {
            "method": "Identity by State (IBS)",
            "description": "Counts shared alleles between user and ancient samples",
            "confidence_explanation": {
                "high": f"{CONFIDENCE_HIGH}+ shared SNPs - Strong statistical power",
                "medium": f"{CONFIDENCE_MEDIUM}-{CONFIDENCE_HIGH-1} shared SNPs - Good for patterns",
                "low": f"{CONFIDENCE_LOW}-{CONFIDENCE_MEDIUM-1} shared SNPs - Interpret with caution",
                "very_low": f"<{CONFIDENCE_LOW} shared SNPs - Limited reliability"
            },
            "limitations": [
                "Consumer arrays test ~0.03% of genome",
                "Ancient DNA has missing/degraded data",
                "High similarity ≠ direct descent",
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
    top_5 = []
    for m in json_data.get("top_matches", [])[:5]:
        top_5.append({
            "name": m["name"],
            "culture": m["culture_name"],
            "period": m["period"],
            "similarity": m["similarity_percent"],
            "shared_snps": m["shared_snps"],
            "confidence": m["confidence"],
            "rank": m.get("rank_display", ""),
            "pmid": m.get("pmid", "")
        })
    
    top_cultures = []
    for c in json_data.get("culture_affinities", [])[:5]:
        top_cultures.append({
            "name": c["name"],
            "affinity": c["affinity_percent"],
            "individuals_compared": c["individuals_compared"],
            "avg_shared_snps": c.get("avg_shared_snps", 0),
            "confidence": c["confidence"],
            "icon": c.get("icon", "🏛️")
        })
    
    stats = json_data.get("statistics", {})
    
    return {
        "top_ancient_matches": top_5,
        "cultural_affinities": top_cultures,
        "total_compared": stats.get("total_individuals_compared", 0),
        "avg_snps_per_comparison": stats.get("avg_snps_compared", 0),
        "high_confidence_matches": stats.get("high_confidence_matches", 0),
        "best_match": {
            "name": stats.get("best_match_name", "N/A"),
            "similarity": stats.get("best_match_similarity", 0)
        },
        "methodology_note": "IBS-based genetic distance. Confidence based on shared SNP count. Educational, not precise ancestry.",
        "full_data": json_data
    }
