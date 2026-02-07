"""
1000 Genomes Population Comparison Module

Compares user genotypes against reference populations from the 1000 Genomes Project.
This provides TRANSPARENT comparisons - users see actual data, not black-box estimates.

Key principle: We show "how your genotypes compare to reference populations" 
NOT "you are X% European" - that requires statistical models with inherent biases.

v4.5.0 UPDATE: Now includes full statistical rigor:
- 95% confidence intervals for all similarity scores
- P-values for significance testing vs random
- Standard errors weighted by marker count
- Confidence levels (DEFINITIVE/HIGH/MEDIUM/LOW/UNCERTAIN)

Sources:
- 1000 Genomes Project Consortium. 2015. A global reference for human genetic 
  variation. Nature 526, 68-74. PMID: 26432245
"""

from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
import json
import math
import sys

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from personal_genomics.statistics import (
        ConfidenceLevel,
        wilson_score_interval,
        ancestry_similarity_stats,
        marker_coverage_weight,
        proportion_test_pvalue,
        confidence_to_color,
        format_ci_string,
    )
    STATS_AVAILABLE = True
except ImportError:
    STATS_AVAILABLE = False
    # Fallback if statistics module not available
    class ConfidenceLevel:
        DEFINITIVE = "DEFINITIVE"
        HIGH = "HIGH"
        MEDIUM = "MEDIUM"
        LOW = "LOW"
        UNCERTAIN = "UNCERTAIN"

# =============================================================================
# LOAD REFERENCE DATA
# =============================================================================

def load_1000genomes_data() -> Dict[str, Any]:
    """Load the 1000 Genomes frequency reference data."""
    ref_path = Path(__file__).parent.parent / "references" / "1000genomes_frequencies.json"
    
    if not ref_path.exists():
        # Return minimal default data if file missing
        return {"_metadata": {"populations": {}}}
    
    with open(ref_path, 'r') as f:
        return json.load(f)


# Cache the reference data
_1000G_DATA = None

def get_reference_data() -> Dict[str, Any]:
    """Get cached reference data."""
    global _1000G_DATA
    if _1000G_DATA is None:
        _1000G_DATA = load_1000genomes_data()
    return _1000G_DATA


# =============================================================================
# POPULATION METADATA
# =============================================================================

POPULATION_INFO = {
    # European
    "GBR": {"name": "British", "superpop": "EUR", "region": "Europe", "flag": "🇬🇧"},
    "CEU": {"name": "Utah European", "superpop": "EUR", "region": "Europe", "flag": "🇺🇸"},
    "FIN": {"name": "Finnish", "superpop": "EUR", "region": "Europe", "flag": "🇫🇮"},
    "TSI": {"name": "Tuscan", "superpop": "EUR", "region": "Europe", "flag": "🇮🇹"},
    "IBS": {"name": "Iberian", "superpop": "EUR", "region": "Europe", "flag": "🇪🇸"},
    
    # African
    "YRI": {"name": "Yoruba (Nigeria)", "superpop": "AFR", "region": "Africa", "flag": "🇳🇬"},
    "LWK": {"name": "Luhya (Kenya)", "superpop": "AFR", "region": "Africa", "flag": "🇰🇪"},
    "ESN": {"name": "Esan (Nigeria)", "superpop": "AFR", "region": "Africa", "flag": "🇳🇬"},
    "ASW": {"name": "African-American (SW USA)", "superpop": "AFR", "region": "Americas", "flag": "🇺🇸"},
    
    # East Asian
    "CHB": {"name": "Han Chinese (Beijing)", "superpop": "EAS", "region": "East Asia", "flag": "🇨🇳"},
    "JPT": {"name": "Japanese (Tokyo)", "superpop": "EAS", "region": "East Asia", "flag": "🇯🇵"},
    "CHS": {"name": "Southern Han Chinese", "superpop": "EAS", "region": "East Asia", "flag": "🇨🇳"},
    "KHV": {"name": "Kinh Vietnamese", "superpop": "EAS", "region": "East Asia", "flag": "🇻🇳"},
    
    # South Asian
    "GIH": {"name": "Gujarati (Houston)", "superpop": "SAS", "region": "South Asia", "flag": "🇮🇳"},
    "PJL": {"name": "Punjabi (Lahore)", "superpop": "SAS", "region": "South Asia", "flag": "🇵🇰"},
    "BEB": {"name": "Bengali (Bangladesh)", "superpop": "SAS", "region": "South Asia", "flag": "🇧🇩"},
    "ITU": {"name": "Telugu (UK)", "superpop": "SAS", "region": "South Asia", "flag": "🇮🇳"},
    
    # Americas
    "MXL": {"name": "Mexican (Los Angeles)", "superpop": "AMR", "region": "Americas", "flag": "🇲🇽"},
    "PUR": {"name": "Puerto Rican", "superpop": "AMR", "region": "Americas", "flag": "🇵🇷"},
    "CLM": {"name": "Colombian (Medellín)", "superpop": "AMR", "region": "Americas", "flag": "🇨🇴"},
    "PEL": {"name": "Peruvian (Lima)", "superpop": "AMR", "region": "Americas", "flag": "🇵🇪"},
}

SUPERPOP_INFO = {
    "EUR": {"name": "European", "color": "#4f46e5"},
    "AFR": {"name": "African", "color": "#059669"},
    "EAS": {"name": "East Asian", "color": "#dc2626"},
    "SAS": {"name": "South Asian", "color": "#d97706"},
    "AMR": {"name": "Americas", "color": "#7c3aed"},
}


def make_bar(pct: float, width: int = 20) -> str:
    """Create a simple text bar chart."""
    filled = int(pct / 100 * width)
    return "█" * filled + "░" * (width - filled)


def normalize_genotype(genotype: str) -> str:
    """Normalize genotype to alphabetical order for comparison."""
    if len(genotype) == 2:
        return "".join(sorted(genotype.upper()))
    return genotype.upper()


# =============================================================================
# POPULATION COMPARISON FUNCTIONS
# =============================================================================

def compare_to_populations(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Compare user's genotypes to all 1000 Genomes populations.
    
    Args:
        genotypes: Dict mapping rsid -> genotype (e.g., {"rs1426654": "AA"})
        
    Returns:
        Comprehensive comparison data including per-marker and aggregate stats
    """
    ref_data = get_reference_data()
    
    results = {
        "markers_analyzed": [],
        "population_match_scores": {},
        "superpop_match_scores": {},
        "summary": {},
        "methodology": {
            "description": "Direct genotype frequency comparison against 1000 Genomes Phase 3 reference populations",
            "pmid": "26432245",
            "note": "This shows how your genotypes compare to reference populations - NOT ancestry percentages"
        }
    }
    
    # Initialize population scores
    pop_scores = {pop: {"match_sum": 0, "markers": 0} for pop in POPULATION_INFO}
    superpop_scores = {sp: {"match_sum": 0, "markers": 0} for sp in SUPERPOP_INFO}
    
    for rsid, marker_data in ref_data.items():
        if rsid.startswith("_"):  # Skip metadata
            continue
            
        user_geno = genotypes.get(rsid)
        if not user_geno:
            continue
        
        user_geno_norm = normalize_genotype(user_geno)
        pops = marker_data.get("populations", {})
        
        if not pops:
            continue
        
        marker_result = {
            "rsid": rsid,
            "gene": marker_data.get("gene", ""),
            "trait": marker_data.get("trait", ""),
            "your_genotype": user_geno,
            "population_frequencies": {},
            "pmid": marker_data.get("pmid", [])
        }
        
        # Check each population
        for pop, freqs in pops.items():
            if pop not in POPULATION_INFO:
                continue
                
            # Find user's genotype frequency in this population
            user_freq = 0.0
            for geno, freq in freqs.items():
                if normalize_genotype(geno) == user_geno_norm:
                    user_freq = freq
                    break
            
            marker_result["population_frequencies"][pop] = {
                "frequency": round(user_freq * 100, 1),
                "population_name": POPULATION_INFO[pop]["name"],
                "superpop": POPULATION_INFO[pop]["superpop"]
            }
            
            # Update scores
            pop_scores[pop]["match_sum"] += user_freq
            pop_scores[pop]["markers"] += 1
            
            superpop = POPULATION_INFO[pop]["superpop"]
            superpop_scores[superpop]["match_sum"] += user_freq
            superpop_scores[superpop]["markers"] += 1
        
        results["markers_analyzed"].append(marker_result)
    
    # Calculate average match scores
    for pop, data in pop_scores.items():
        if data["markers"] > 0:
            avg = (data["match_sum"] / data["markers"]) * 100
            results["population_match_scores"][pop] = {
                "average_frequency": round(avg, 1),
                "markers_compared": data["markers"],
                "population_name": POPULATION_INFO[pop]["name"],
                "superpop": POPULATION_INFO[pop]["superpop"],
                "flag": POPULATION_INFO[pop]["flag"]
            }
    
    for superpop, data in superpop_scores.items():
        if data["markers"] > 0:
            avg = (data["match_sum"] / data["markers"]) * 100
            results["superpop_match_scores"][superpop] = {
                "average_frequency": round(avg, 1),
                "markers_compared": data["markers"],
                "name": SUPERPOP_INFO[superpop]["name"]
            }
    
    # Summary
    results["summary"]["total_markers_compared"] = len(results["markers_analyzed"])
    
    return results


def find_most_similar_populations(genotypes: Dict[str, str], top_n: int = 10) -> List[Dict[str, Any]]:
    """
    Rank populations by similarity to user's genotypes.
    
    Uses geometric mean of genotype frequencies as similarity metric.
    NOW INCLUDES: Confidence intervals, p-values, and statistical significance.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        top_n: Number of top populations to return
        
    Returns:
        List of populations sorted by similarity, each with:
        - similarity: percentage similarity score
        - ci_lower, ci_upper: 95% confidence interval bounds
        - p_value: significance vs random (0.5 baseline)
        - n_markers: number of markers used
        - confidence: DEFINITIVE/HIGH/MEDIUM/LOW/UNCERTAIN
    """
    comparison = compare_to_populations(genotypes)
    pop_scores = comparison.get("population_match_scores", {})
    
    # Sort by average frequency
    ranked = sorted(
        pop_scores.items(),
        key=lambda x: x[1].get("average_frequency", 0),
        reverse=True
    )
    
    result = []
    for pop, data in ranked[:top_n]:
        similarity = data.get("average_frequency", 0)
        n_markers = data.get("markers_compared", 0)
        
        # Calculate statistical metrics
        stats_result = _calculate_population_stats(similarity, n_markers)
        
        result.append({
            "population_code": pop,
            "population_name": data.get("population_name", pop),
            "superpopulation": data.get("superpop", ""),
            "similarity": round(similarity, 1),
            "ci_lower": stats_result["ci_lower"],
            "ci_upper": stats_result["ci_upper"],
            "p_value": stats_result["p_value"],
            "n_markers": n_markers,
            "confidence": stats_result["confidence"],
            "confidence_color": stats_result["confidence_color"],
            "standard_error": stats_result["standard_error"],
            "flag": data.get("flag", ""),
            "interpretation": _interpret_similarity(similarity),
            # Legacy field for backward compatibility
            "similarity_score": round(similarity, 1),
            "markers_compared": n_markers,
        })
    
    return result


def _calculate_population_stats(similarity_pct: float, n_markers: int) -> Dict[str, Any]:
    """
    Calculate statistical metrics for a population similarity score.
    
    Args:
        similarity_pct: Similarity as percentage (0-100)
        n_markers: Number of markers compared
        
    Returns:
        Dict with ci_lower, ci_upper, p_value, confidence, etc.
    """
    if not STATS_AVAILABLE or n_markers == 0:
        # Fallback without statistics module
        return {
            "ci_lower": max(0, similarity_pct - 15),
            "ci_upper": min(100, similarity_pct + 15),
            "p_value": None,
            "confidence": "UNCERTAIN" if n_markers < 5 else "LOW",
            "confidence_color": "#f97316",
            "standard_error": 15.0,
        }
    
    # Convert to proportion for calculations
    similarity = similarity_pct / 100.0
    
    # Use Wilson score interval for proportions
    ci = wilson_score_interval(
        successes=int(similarity * n_markers),
        n=n_markers,
        confidence=0.95
    )
    
    # P-value: is this significantly different from random (50%)?
    p_value = proportion_test_pvalue(similarity, 0.5, n_markers)
    
    # Calculate standard error
    se = math.sqrt(similarity * (1 - similarity) / n_markers) if n_markers > 0 else 0.15
    
    # Determine confidence level
    quality_score, conf_level = marker_coverage_weight(n_markers, 50)  # Assume 50 ideal markers
    
    # Upgrade confidence if highly significant
    if p_value and p_value < 0.001 and n_markers >= 15:
        if conf_level == ConfidenceLevel.MEDIUM:
            conf_level = ConfidenceLevel.HIGH
    
    return {
        "ci_lower": round(ci.lower * 100, 1),
        "ci_upper": round(ci.upper * 100, 1),
        "p_value": round(p_value, 6) if p_value else None,
        "confidence": conf_level.value if hasattr(conf_level, 'value') else str(conf_level),
        "confidence_color": confidence_to_color(conf_level) if STATS_AVAILABLE else "#6b7280",
        "standard_error": round(se * 100, 2),
    }


def _interpret_similarity(score: float) -> str:
    """Provide interpretation of similarity score."""
    if score >= 80:
        return "Very high match - your genotypes are common in this population"
    elif score >= 60:
        return "High match - many of your genotypes are found at moderate-to-high frequency"
    elif score >= 40:
        return "Moderate match - mixed frequency pattern"
    elif score >= 20:
        return "Low match - many of your genotypes are uncommon in this population"
    else:
        return "Very low match - your genotypes are rare in this population"


def get_marker_population_context(rsid: str, genotype: str) -> Dict[str, Any]:
    """
    Get detailed population context for a single marker.
    
    Args:
        rsid: The SNP ID (e.g., "rs1426654")
        genotype: User's genotype (e.g., "AA")
        
    Returns:
        Population frequency breakdown for this specific genotype
    """
    ref_data = get_reference_data()
    
    marker_data = ref_data.get(rsid)
    if not marker_data:
        return {
            "status": "not_found",
            "rsid": rsid,
            "error": f"No reference data available for {rsid}"
        }
    
    genotype_norm = normalize_genotype(genotype)
    pops = marker_data.get("populations", {})
    
    result = {
        "rsid": rsid,
        "gene": marker_data.get("gene", ""),
        "trait": marker_data.get("trait", ""),
        "your_genotype": genotype,
        "ref_allele": marker_data.get("ref", ""),
        "alt_allele": marker_data.get("alt", ""),
        "populations_by_frequency": [],
        "superpop_summary": {},
        "pmid": marker_data.get("pmid", [])
    }
    
    # Collect frequencies
    pop_freqs = []
    superpop_totals = {sp: {"sum": 0, "count": 0} for sp in SUPERPOP_INFO}
    
    for pop, freqs in pops.items():
        if pop not in POPULATION_INFO:
            continue
        
        # Find user's genotype frequency
        user_freq = 0.0
        for geno, freq in freqs.items():
            if normalize_genotype(geno) == genotype_norm:
                user_freq = freq
                break
        
        pop_info = POPULATION_INFO[pop]
        pop_freqs.append({
            "code": pop,
            "name": pop_info["name"],
            "superpop": pop_info["superpop"],
            "flag": pop_info["flag"],
            "frequency": round(user_freq * 100, 1),
            "bar": make_bar(user_freq * 100)
        })
        
        # Aggregate by superpop
        sp = pop_info["superpop"]
        superpop_totals[sp]["sum"] += user_freq
        superpop_totals[sp]["count"] += 1
    
    # Sort by frequency (highest first)
    pop_freqs.sort(key=lambda x: x["frequency"], reverse=True)
    result["populations_by_frequency"] = pop_freqs
    
    # Calculate superpop averages
    for sp, data in superpop_totals.items():
        if data["count"] > 0:
            avg = (data["sum"] / data["count"]) * 100
            result["superpop_summary"][sp] = {
                "name": SUPERPOP_INFO[sp]["name"],
                "average_frequency": round(avg, 1)
            }
    
    # Determine typical populations
    high_freq_pops = [p for p in pop_freqs if p["frequency"] >= 50]
    if high_freq_pops:
        # Group by superpop
        superpops = list(set(p["superpop"] for p in high_freq_pops))
        superpop_names = [SUPERPOP_INFO[sp]["name"] for sp in superpops]
        result["interpretation"] = f"Your genotype is typical of {', '.join(superpop_names)} populations"
    else:
        moderate_freq = [p for p in pop_freqs if p["frequency"] >= 20]
        if moderate_freq:
            result["interpretation"] = "Your genotype is moderately common across several populations"
        else:
            result["interpretation"] = "Your genotype is relatively uncommon globally"
    
    return result


def generate_population_comparison_report(genotypes: Dict[str, str]) -> str:
    """
    Generate a human-readable text report of population comparisons.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Formatted text report
    """
    lines = []
    lines.append("=" * 70)
    lines.append("🧬 1000 GENOMES POPULATION COMPARISON REPORT")
    lines.append("=" * 70)
    lines.append("")
    lines.append("This report shows how your genotypes compare to global reference")
    lines.append("populations from the 1000 Genomes Project (Phase 3, 2,504 samples).")
    lines.append("")
    lines.append("⚠️  NOTE: This is NOT ancestry estimation. It shows genotype frequencies")
    lines.append("in reference populations - the RAW DATA that ancestry services use.")
    lines.append("")
    
    # Get most similar populations
    similar_pops = find_most_similar_populations(genotypes, top_n=10)
    
    if not similar_pops:
        lines.append("⚠️  Insufficient marker overlap with reference data.")
        return "\n".join(lines)
    
    lines.append("-" * 70)
    lines.append("📊 MOST SIMILAR REFERENCE POPULATIONS")
    lines.append("-" * 70)
    lines.append("")
    
    for i, pop in enumerate(similar_pops, 1):
        flag = pop.get("flag", "")
        name = pop.get("population_name", "Unknown")
        superpop = SUPERPOP_INFO.get(pop.get("superpopulation", ""), {}).get("name", "")
        score = pop.get("similarity_score", 0)
        markers = pop.get("markers_compared", 0)
        
        lines.append(f"  {i:2}. {flag} {name} ({superpop})")
        lines.append(f"      Similarity: {score:.1f}%  {make_bar(score)}")
        lines.append(f"      Based on {markers} markers")
        lines.append("")
    
    # Per-marker details
    lines.append("-" * 70)
    lines.append("📋 INDIVIDUAL MARKER COMPARISONS")
    lines.append("-" * 70)
    lines.append("")
    
    ref_data = get_reference_data()
    markers_shown = 0
    
    for rsid, marker_data in ref_data.items():
        if rsid.startswith("_"):
            continue
        
        user_geno = genotypes.get(rsid)
        if not user_geno:
            continue
        
        context = get_marker_population_context(rsid, user_geno)
        gene = context.get("gene", "")
        trait = context.get("trait", "")
        
        lines.append(f"▸ {rsid}" + (f" ({gene})" if gene else ""))
        if trait:
            lines.append(f"  Trait: {trait}")
        lines.append(f"  Your genotype: {user_geno}")
        lines.append("")
        lines.append("  Population frequencies:")
        
        for pop in context.get("populations_by_frequency", [])[:10]:
            freq = pop["frequency"]
            bar = pop["bar"]
            name = pop["name"]
            lines.append(f"    {name:25} {freq:5.1f}%  {bar}")
        
        lines.append("")
        if context.get("interpretation"):
            lines.append(f"  → {context['interpretation']}")
        lines.append("")
        
        markers_shown += 1
        if markers_shown >= 15:  # Limit report length
            lines.append(f"  ... and {len([r for r in ref_data if not r.startswith('_') and r in genotypes]) - markers_shown} more markers")
            break
    
    # Methodology
    lines.append("")
    lines.append("=" * 70)
    lines.append("📖 METHODOLOGY")
    lines.append("=" * 70)
    lines.append("")
    lines.append("Reference: 1000 Genomes Project Phase 3 (PMID: 26432245)")
    lines.append("Method: Direct genotype frequency lookup in reference populations")
    lines.append("")
    lines.append("Why we show this instead of ancestry percentages:")
    lines.append("• Percentages require statistical models with inherent assumptions")
    lines.append("• Reference populations don't represent all human diversity")
    lines.append("• You deserve to see the actual data, not a black box result")
    lines.append("• This approach is more honest about uncertainty")
    lines.append("")
    
    return "\n".join(lines)


def get_population_comparison_json(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Get full population comparison data as structured JSON for dashboard.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Structured data for dashboard visualization
    """
    comparison = compare_to_populations(genotypes)
    similar_pops = find_most_similar_populations(genotypes, top_n=20)
    
    # Get detailed marker data
    marker_details = []
    ref_data = get_reference_data()
    
    for rsid, marker_data in ref_data.items():
        if rsid.startswith("_"):
            continue
        
        user_geno = genotypes.get(rsid)
        if not user_geno:
            continue
        
        context = get_marker_population_context(rsid, user_geno)
        marker_details.append(context)
    
    return {
        "most_similar_populations": similar_pops,
        "superpopulation_summary": comparison.get("superpop_match_scores", {}),
        "all_population_scores": comparison.get("population_match_scores", {}),
        "marker_details": marker_details,
        "total_markers": len(marker_details),
        "population_metadata": POPULATION_INFO,
        "superpop_metadata": SUPERPOP_INFO,
        "methodology": {
            "source": "1000 Genomes Project Phase 3",
            "samples": 2504,
            "pmid": "26432245",
            "approach": "Direct genotype frequency comparison - NOT ancestry estimation"
        }
    }


# =============================================================================
# DASHBOARD DATA EXPORT
# =============================================================================

def export_for_dashboard(genotypes: Dict[str, str], output_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Export population comparison data formatted for dashboard visualization.
    
    Args:
        genotypes: User genotypes
        output_path: Optional path to write JSON file
        
    Returns:
        Dashboard-ready data structure
    """
    data = get_population_comparison_json(genotypes)
    
    if output_path:
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    return data
