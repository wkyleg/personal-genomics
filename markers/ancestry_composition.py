"""
Ancestry Composition Analysis Module
Population admixture estimation from SNP data

This module provides:
- Reference population comparisons at CONTINENTAL level only
- Admixture estimation using ancestry informative markers (AIMs)
- Confidence intervals for all estimates
- Explicit methodology limitations

╔══════════════════════════════════════════════════════════════════════════════╗
║  ⚠️  IMPORTANT LIMITATIONS                                                   ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Consumer DNA arrays have significant limitations for ancestry estimation:   ║
║                                                                              ║
║  • Can only reliably distinguish CONTINENTAL ancestry (not sub-regions)      ║
║  • Cannot distinguish Irish from Scottish, German from Polish, etc.          ║
║  • Estimates have wide confidence intervals (typically ±10-20%)              ║
║  • Reference panels are biased toward well-studied populations               ║
║  • Ancient ancestry may not match modern population labels                   ║
║                                                                              ║
║  These results provide a rough estimate only. For genealogical research,     ║
║  combine with paper trail documentation.                                     ║
╚══════════════════════════════════════════════════════════════════════════════╝

Sources & PMIDs:
- 1000 Genomes Project Consortium. 2015. A global reference for human genetic 
  variation. PMID: 26432245
- The International HapMap 3 Consortium. 2010. Integrating common and rare 
  genetic variation. PMID: 20981092
- Rosenberg NA et al. 2002. Genetic structure of human populations. 
  PMID: 12493913
- Li JZ et al. 2008. Worldwide human relationships inferred from genome-wide 
  patterns. PMID: 18292342
"""

from typing import Dict, List, Optional, Any, Tuple
import math

# =============================================================================
# METHODOLOGY DISCLAIMER
# =============================================================================

ANCESTRY_DISCLAIMER = """
⚠️ ANCESTRY ESTIMATION LIMITATIONS

This ancestry analysis can only reliably estimate CONTINENTAL-level ancestry:
• European (EUR)
• African (AFR)
• East Asian (EAS)
• South Asian (SAS)
• Americas/Indigenous (AMR)

SUB-REGIONAL claims (Irish vs Scottish, Nigerian vs Ghanaian) are NOT reliable
with consumer DNA arrays due to:
1. Insufficient marker density for fine-scale population structure
2. Reference panel limitations
3. Recent shared ancestry between neighboring populations

All percentages include confidence intervals showing the range of uncertainty.
"""

# =============================================================================
# ANCESTRY INFORMATIVE MARKERS (AIMs)
# =============================================================================

# Population frequency data from 1000 Genomes (PMID: 26432245) and published studies
# Format: rsid -> {populations: {pop: allele_freq}, ancestry_info, pmid}

ANCESTRY_INFORMATIVE_MARKERS = {
    # =========================================================================
    # CONTINENTAL-LEVEL MARKERS (High differentiation)
    # =========================================================================
    
    # SLC24A5 - Major skin pigmentation gene (European selection)
    "rs1426654": {
        "gene": "SLC24A5",
        "trait": "skin_pigmentation",
        "ancestral": "G",
        "derived": "A",
        "frequencies": {
            "EUR": 0.98,
            "AFR": 0.03,
            "EAS": 0.03,
            "SAS": 0.85,
            "AMR": 0.55,
        },
        "delta": 0.95,
        "informativeness": "continental",
        "pmid": ["16357253", "26432245"]
    },
    
    # SLC45A2 - Skin/hair pigmentation
    "rs16891982": {
        "gene": "SLC45A2",
        "trait": "skin_pigmentation",
        "ancestral": "C",
        "derived": "G",
        "frequencies": {
            "EUR": 0.96,
            "AFR": 0.02,
            "EAS": 0.02,
            "SAS": 0.15,
            "AMR": 0.45,
        },
        "delta": 0.94,
        "informativeness": "continental",
        "pmid": ["17999355", "26432245"]
    },
    
    # HERC2/OCA2 - Eye color (European)
    "rs12913832": {
        "gene": "HERC2",
        "trait": "eye_color",
        "ancestral": "A",
        "derived": "G",
        "frequencies": {
            "EUR": 0.75,
            "AFR": 0.01,
            "EAS": 0.01,
            "SAS": 0.10,
            "AMR": 0.35,
        },
        "delta": 0.74,
        "informativeness": "european_indicator",
        "pmid": ["18252222", "26432245"]
    },
    
    # LCT - Lactase persistence (European, some African)
    "rs4988235": {
        "gene": "LCT",
        "trait": "lactase_persistence",
        "ancestral": "G",
        "derived": "A",
        "frequencies": {
            "EUR": 0.75,
            "AFR": 0.10,
            "EAS": 0.02,
            "SAS": 0.25,
            "AMR": 0.40,
        },
        "delta": 0.73,
        "informativeness": "european_indicator",
        "pmid": ["11788828", "26432245"]
    },
    
    # EDAR - Hair morphology (East Asian)
    "rs3827760": {
        "gene": "EDAR",
        "trait": "hair_morphology",
        "ancestral": "A",
        "derived": "G",
        "frequencies": {
            "EUR": 0.01,
            "AFR": 0.00,
            "EAS": 0.93,
            "SAS": 0.05,
            "AMR": 0.65,
        },
        "delta": 0.93,
        "informativeness": "east_asian_specific",
        "pmid": ["18561325", "26432245"]
    },
    
    # ABCC11 - Earwax/body odor (East Asian)
    "rs17822931": {
        "gene": "ABCC11",
        "trait": "earwax_type",
        "ancestral": "G",
        "derived": "A",
        "frequencies": {
            "EUR": 0.12,
            "AFR": 0.02,
            "EAS": 0.95,
            "SAS": 0.15,
            "AMR": 0.55,
        },
        "delta": 0.93,
        "informativeness": "east_asian_specific",
        "pmid": ["16444273", "26432245"]
    },
    
    # Duffy antigen - African malaria adaptation
    "rs2814778": {
        "gene": "DARC",
        "trait": "duffy_null",
        "ancestral": "T",
        "derived": "C",
        "frequencies": {
            "EUR": 0.00,
            "AFR": 0.97,
            "EAS": 0.00,
            "SAS": 0.00,
            "AMR": 0.15,
        },
        "delta": 0.97,
        "informativeness": "african_specific",
        "pmid": ["11171069", "26432245"]
    },
    
    # =========================================================================
    # ADDITIONAL CONTINENTAL MARKERS
    # =========================================================================
    
    "rs1800414": {
        "gene": "OCA2",
        "trait": "pigmentation",
        "frequencies": {
            "EUR": 0.05,
            "AFR": 0.00,
            "EAS": 0.65,
            "SAS": 0.10,
            "AMR": 0.35,
        },
        "delta": 0.65,
        "informativeness": "east_asian_indicator",
        "pmid": ["17952075", "26432245"]
    },
    
    "rs2228479": {
        "gene": "MC1R",
        "trait": "pigmentation",
        "frequencies": {
            "EUR": 0.12,
            "AFR": 0.00,
            "EAS": 0.00,
            "SAS": 0.02,
            "AMR": 0.05,
        },
        "delta": 0.12,
        "informativeness": "european_indicator",
        "pmid": ["11017081", "26432245"]
    },
    
    "rs885479": {
        "gene": "MC1R",
        "trait": "pigmentation",
        "frequencies": {
            "EUR": 0.03,
            "AFR": 0.00,
            "EAS": 0.65,
            "SAS": 0.15,
            "AMR": 0.30,
        },
        "delta": 0.65,
        "informativeness": "east_asian_indicator",
        "pmid": ["11017081", "26432245"]
    },
    
    # Native American informative
    "rs3811801": {
        "gene": "FUT2",
        "trait": "secretor_status",
        "frequencies": {
            "EUR": 0.45,
            "AFR": 0.35,
            "EAS": 0.40,
            "SAS": 0.40,
            "AMR": 0.10,
        },
        "delta": 0.35,
        "informativeness": "native_american_indicator",
        "pmid": ["26432245"]
    },
    
    # South Asian informative markers
    "rs1042602": {
        "gene": "TYR",
        "trait": "pigmentation",
        "frequencies": {
            "EUR": 0.40,
            "AFR": 0.05,
            "EAS": 0.02,
            "SAS": 0.30,
            "AMR": 0.25,
        },
        "delta": 0.38,
        "informativeness": "subcontinental",
        "pmid": ["18483556", "26432245"]
    },
    
    "rs260690": {
        "gene": "TYRP1",
        "frequencies": {"EUR": 0.50, "AFR": 0.08, "EAS": 0.15, "SAS": 0.25, "AMR": 0.30},
        "delta": 0.42,
        "informativeness": "continental",
        "pmid": ["26432245"]
    },
    
    "rs1408799": {
        "gene": "TYRP1",
        "frequencies": {"EUR": 0.75, "AFR": 0.15, "EAS": 0.20, "SAS": 0.35, "AMR": 0.45},
        "delta": 0.60,
        "informativeness": "continental",
        "pmid": ["26432245"]
    },
    
    "rs7495174": {
        "gene": "OCA2",
        "frequencies": {"EUR": 0.85, "AFR": 0.02, "EAS": 0.65, "SAS": 0.50, "AMR": 0.55},
        "delta": 0.83,
        "informativeness": "continental",
        "pmid": ["26432245"]
    },
    
    "rs4778138": {
        "gene": "OCA2",
        "frequencies": {"EUR": 0.80, "AFR": 0.15, "EAS": 0.65, "SAS": 0.45, "AMR": 0.50},
        "delta": 0.65,
        "informativeness": "continental",
        "pmid": ["26432245"]
    },
    
    "rs1393350": {
        "gene": "TYR",
        "frequencies": {"EUR": 0.25, "AFR": 0.05, "EAS": 0.01, "SAS": 0.08, "AMR": 0.12},
        "delta": 0.24,
        "informativeness": "european_indicator",
        "pmid": ["26432245"]
    },
    
    "rs2305498": {
        "frequencies": {"EUR": 0.30, "AFR": 0.65, "EAS": 0.25, "SAS": 0.35, "AMR": 0.40},
        "delta": 0.40,
        "informativeness": "african_indicator",
        "pmid": ["26432245"]
    },
    
    "rs10843104": {
        "frequencies": {"EUR": 0.15, "AFR": 0.75, "EAS": 0.10, "SAS": 0.20, "AMR": 0.35},
        "delta": 0.65,
        "informativeness": "african_indicator",
        "pmid": ["26432245"]
    },
    
    "rs2789823": {
        "frequencies": {"EUR": 0.45, "AFR": 0.10, "EAS": 0.85, "SAS": 0.30, "AMR": 0.50},
        "delta": 0.75,
        "informativeness": "east_asian_indicator",
        "pmid": ["26432245"]
    },
    
    "rs2070959": {
        "gene": "UGT1A1",
        "frequencies": {"EUR": 0.35, "AFR": 0.45, "EAS": 0.12, "SAS": 0.25, "AMR": 0.30},
        "delta": 0.33,
        "informativeness": "subcontinental",
        "pmid": ["26432245"]
    },
}

# =============================================================================
# POPULATION REFERENCE DATA
# =============================================================================

POPULATION_DESCRIPTIONS = {
    "EUR": {
        "name": "European",
        "abbreviation": "EUR",
        "description": "European ancestry - cannot distinguish sub-regions (Northern, Southern, Eastern, Western)",
        "limitations": [
            "Cannot distinguish British from German from Italian",
            "Cannot identify specific European country of origin",
            "Reflects genetic similarity, not cultural/national identity"
        ],
        "typical_haplogroups_mtdna": ["H", "U", "K", "J", "T", "V"],
        "typical_haplogroups_ydna": ["R1b", "R1a", "I", "J2", "E1b1b"],
        "pmid": ["26432245", "12493913"]
    },
    "AFR": {
        "name": "African",
        "abbreviation": "AFR",
        "description": "Sub-Saharan African ancestry - cannot distinguish sub-regions",
        "limitations": [
            "Cannot distinguish West African from East African",
            "Cannot identify specific ethnic groups or countries",
            "African genetic diversity is highest globally"
        ],
        "typical_haplogroups_mtdna": ["L0", "L1", "L2", "L3"],
        "typical_haplogroups_ydna": ["E", "A", "B"],
        "pmid": ["26432245", "12493913"]
    },
    "EAS": {
        "name": "East Asian",
        "abbreviation": "EAS",
        "description": "East Asian ancestry - cannot distinguish Chinese, Japanese, Korean, etc.",
        "limitations": [
            "Cannot distinguish Chinese from Japanese from Korean",
            "Cannot identify specific sub-populations",
            "Southeast Asian may appear as partial East Asian"
        ],
        "typical_haplogroups_mtdna": ["D", "B", "F", "A", "C", "M"],
        "typical_haplogroups_ydna": ["O", "C", "D", "N"],
        "pmid": ["26432245", "12493913"]
    },
    "SAS": {
        "name": "South Asian",
        "abbreviation": "SAS",
        "description": "South Asian ancestry - cannot distinguish specific regions",
        "limitations": [
            "Cannot distinguish Indian from Pakistani from Bangladeshi",
            "Cannot identify caste or ethnic group",
            "High internal diversity in South Asia"
        ],
        "typical_haplogroups_mtdna": ["M", "R", "U"],
        "typical_haplogroups_ydna": ["R1a", "H", "L", "J2"],
        "pmid": ["26432245", "12493913"]
    },
    "AMR": {
        "name": "Americas (Indigenous + Admixed)",
        "abbreviation": "AMR",
        "description": "Indigenous American ancestry component - highly variable",
        "limitations": [
            "Most individuals with AMR ancestry are admixed",
            "Cannot distinguish specific Indigenous groups",
            "Pre-Columbian diversity not well captured"
        ],
        "typical_haplogroups_mtdna": ["A", "B", "C", "D", "X"],
        "typical_haplogroups_ydna": ["Q", "C"],
        "pmid": ["26432245", "12493913"]
    },
}


def calculate_wilson_confidence_interval(
    successes: int,
    trials: int,
    confidence: float = 0.95
) -> Tuple[float, float]:
    """
    Calculate Wilson score confidence interval for a proportion.
    
    More accurate than normal approximation, especially for extreme proportions.
    
    Args:
        successes: Number of "successes" (matching alleles)
        trials: Total number of trials
        confidence: Confidence level (default 95%)
        
    Returns:
        Tuple of (lower_bound, upper_bound) as proportions
    """
    if trials == 0:
        return (0.0, 1.0)
    
    # Z-score for confidence level
    z = 1.96 if confidence == 0.95 else 1.645 if confidence == 0.90 else 2.576
    
    p_hat = successes / trials
    
    denominator = 1 + z**2 / trials
    center = (p_hat + z**2 / (2 * trials)) / denominator
    
    margin = z * math.sqrt((p_hat * (1 - p_hat) + z**2 / (4 * trials)) / trials) / denominator
    
    lower = max(0, center - margin)
    upper = min(1, center + margin)
    
    return (lower, upper)


def estimate_ancestry(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Estimate ancestry composition from available AIMs.
    
    Returns CONTINENTAL-LEVEL estimates only with confidence intervals.
    Sub-regional ancestry (e.g., Irish vs Scottish) CANNOT be reliably determined.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Dict with ancestry estimates, confidence intervals, and methodology notes
    """
    # Initialize population scores
    populations = ["EUR", "AFR", "EAS", "SAS", "AMR"]
    log_likelihoods = {pop: 0.0 for pop in populations}
    markers_used = 0
    marker_details = []
    all_pmids = set()
    
    for rsid, info in ANCESTRY_INFORMATIVE_MARKERS.items():
        geno = genotypes.get(rsid)
        if not geno:
            continue
        
        freqs = info.get("frequencies", {})
        if not all(pop in freqs for pop in populations):
            continue
        
        markers_used += 1
        
        # Collect PMIDs
        if "pmid" in info:
            all_pmids.update(info["pmid"])
        
        # Determine allele count
        derived = info.get("derived", info.get("ancestral", ""))
        
        derived_count = sum(1 for a in geno if a.upper() == derived.upper())
        
        # Calculate genotype probability for each population
        for pop in populations:
            p = freqs.get(pop, 0.5)
            q = 1 - p
            
            # Avoid log(0)
            p = max(0.001, min(0.999, p))
            q = max(0.001, min(0.999, q))
            
            if derived_count == 2:
                prob = p * p
            elif derived_count == 1:
                prob = 2 * p * q
            else:
                prob = q * q
            
            log_likelihoods[pop] += math.log(prob)
        
        # Track informative markers
        if info.get("delta", 0) > 0.3:
            marker_details.append({
                "rsid": rsid,
                "gene": info.get("gene", ""),
                "genotype": geno,
                "informativeness": info.get("informativeness"),
                "pmid": info.get("pmid", [])
            })
    
    # Convert log-likelihoods to proportions
    if markers_used < 3:
        return {
            "status": "insufficient_data",
            "markers_found": markers_used,
            "error": "Need at least 3 ancestry informative markers for estimation",
            "disclaimer": ANCESTRY_DISCLAIMER
        }
    
    # Normalize to probabilities (softmax-like)
    max_ll = max(log_likelihoods.values())
    exp_lls = {pop: math.exp(ll - max_ll) for pop, ll in log_likelihoods.items()}
    total = sum(exp_lls.values())
    
    # Calculate point estimates (whole numbers only - no false precision)
    proportions = {pop: round(exp_ll / total * 100) for pop, exp_ll in exp_lls.items()}
    
    # Ensure proportions sum to 100
    diff = 100 - sum(proportions.values())
    if diff != 0:
        max_pop = max(proportions.keys(), key=lambda p: proportions[p])
        proportions[max_pop] += diff
    
    # Calculate confidence intervals based on marker count
    # More markers = narrower intervals
    confidence_intervals = {}
    for pop, pct in proportions.items():
        # Use Wilson interval treating proportion as binomial
        # Approximate "successes" based on proportion and marker count
        successes = int(pct / 100 * markers_used * 2)  # diploid
        trials = markers_used * 2
        
        lower, upper = calculate_wilson_confidence_interval(successes, trials)
        
        # Convert to percentage and widen based on low marker count
        # Minimum ±10% uncertainty for ancestry estimates
        interval_width = max(10, int((upper - lower) * 100))
        
        ci_lower = max(0, pct - interval_width // 2)
        ci_upper = min(100, pct + interval_width // 2)
        
        confidence_intervals[pop] = {
            "point_estimate": pct,
            "lower_bound": ci_lower,
            "upper_bound": ci_upper,
            "display": f"{pct}% ({ci_lower}-{ci_upper}%)"
        }
    
    # Determine confidence level
    if markers_used >= 15:
        confidence = "moderate"
        confidence_note = "Sufficient markers for rough continental estimate"
    elif markers_used >= 8:
        confidence = "low"
        confidence_note = "Limited markers - wide confidence intervals"
    else:
        confidence = "very low"
        confidence_note = "Few markers found - estimates highly uncertain"
    
    # Sort by proportion
    sorted_ancestry = sorted(proportions.items(), key=lambda x: x[1], reverse=True)
    
    # Build result
    result = {
        "status": "success",
        "disclaimer": ANCESTRY_DISCLAIMER,
        "markers_used": markers_used,
        "confidence": confidence,
        "confidence_note": confidence_note,
        
        # Ancestry with confidence intervals (continental level only)
        "ancestry_proportions": {
            POPULATION_DESCRIPTIONS[pop]["name"]: confidence_intervals[pop]["display"]
            for pop, pct in sorted_ancestry
            if pct > 0
        },
        
        # Raw data for programmatic use
        "ancestry_proportions_raw": proportions,
        "confidence_intervals": confidence_intervals,
        
        "primary_ancestry": POPULATION_DESCRIPTIONS[sorted_ancestry[0][0]]["name"],
        "primary_percentage": sorted_ancestry[0][1],
        "primary_range": f"{confidence_intervals[sorted_ancestry[0][0]]['lower_bound']}-{confidence_intervals[sorted_ancestry[0][0]]['upper_bound']}%",
        
        "detailed_breakdown": [],
        
        "methodology": {
            "description": "Likelihood-based ancestry estimation using ancestry informative markers (AIMs)",
            "resolution": "CONTINENTAL ONLY - sub-regional ancestry cannot be determined",
            "reference_populations": "1000 Genomes Project superpopulations",
            "marker_count": markers_used,
            "pmid": sorted(list(all_pmids))
        },
        
        "limitations": [
            "Can only distinguish continental ancestry (European, African, East Asian, South Asian, Americas)",
            "CANNOT reliably distinguish sub-regions (e.g., Irish vs Scottish, Nigerian vs Ghanaian)",
            "Confidence intervals are wide - interpret with caution",
            "Reference panels may not represent all populations equally",
            "Modern population labels may not reflect ancient ancestry",
            "Results should be combined with genealogical paper trail research"
        ]
    }
    
    # Add detailed breakdown for significant ancestries (>0%)
    for pop, pct in sorted_ancestry:
        if pct > 0:
            pop_info = POPULATION_DESCRIPTIONS.get(pop, {})
            ci = confidence_intervals[pop]
            result["detailed_breakdown"].append({
                "population": pop_info.get("name", pop),
                "abbreviation": pop,
                "percentage": pct,
                "range": f"{ci['lower_bound']}-{ci['upper_bound']}%",
                "confidence_interval": ci,
                "description": pop_info.get("description", ""),
                "limitations": pop_info.get("limitations", []),
            })
    
    # Key informative markers with PMIDs
    result["informative_markers"] = marker_details[:10]
    
    return result


def detect_admixture(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Detect signs of recent admixture between populations.
    
    Recent admixture shows characteristic patterns:
    - High heterozygosity at AIMs
    - Mixed ancestry percentages in 20-80% range
    """
    ancestry = estimate_ancestry(genotypes)
    
    if ancestry.get("status") != "success":
        return {
            "status": "insufficient_data",
            "is_admixed": "unknown",
            "disclaimer": ANCESTRY_DISCLAIMER
        }
    
    proportions = ancestry.get("ancestry_proportions_raw", {})
    
    # Count populations with significant contribution
    significant_pops = sum(1 for pct in proportions.values() if pct >= 15)
    
    # Check heterozygosity at high-delta AIMs
    het_count = 0
    aim_count = 0
    
    for rsid, info in ANCESTRY_INFORMATIVE_MARKERS.items():
        if info.get("delta", 0) < 0.5:
            continue
            
        geno = genotypes.get(rsid)
        if not geno or len(geno) < 2:
            continue
            
        aim_count += 1
        if geno[0] != geno[1]:
            het_count += 1
    
    het_rate = het_count / aim_count if aim_count > 0 else 0
    
    is_admixed = significant_pops >= 2 or het_rate > 0.6
    
    return {
        "status": "success",
        "disclaimer": ANCESTRY_DISCLAIMER,
        "is_admixed": is_admixed,
        "admixture_confidence": "low" if aim_count < 10 else "moderate",
        "admixture_evidence": {
            "populations_with_15pct_plus": significant_pops,
            "heterozygosity_rate": round(het_rate, 2),
            "aims_checked": aim_count
        },
        "interpretation": _interpret_admixture(proportions, significant_pops, het_rate),
        "ancestry_summary": ancestry.get("ancestry_proportions", {}),
        "methodology_note": "Admixture detection based on AIM heterozygosity and multi-population signal"
    }


def _interpret_admixture(proportions: Dict, sig_pops: int, het_rate: float) -> str:
    """Generate human-readable admixture interpretation."""
    if sig_pops == 1:
        top_pop = max(proportions.items(), key=lambda x: x[1])
        pop_name = POPULATION_DESCRIPTIONS.get(top_pop[0], {}).get('name', top_pop[0])
        return f"Ancestry appears primarily {pop_name} (note: sub-regional origin cannot be determined)"
    elif sig_pops == 2:
        top_two = sorted(proportions.items(), key=lambda x: x[1], reverse=True)[:2]
        pop1 = POPULATION_DESCRIPTIONS.get(top_two[0][0], {}).get('name', top_two[0][0])
        pop2 = POPULATION_DESCRIPTIONS.get(top_two[1][0], {}).get('name', top_two[1][0])
        return f"Possible admixture between {pop1} and {pop2} ancestry (confidence intervals overlap significantly)"
    else:
        return "Complex admixture pattern - multiple continental ancestries detected (wide confidence intervals)"


def get_ancestry_summary(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Complete ancestry analysis with composition and admixture detection.
    
    Returns CONTINENTAL-LEVEL estimates only with confidence intervals.
    """
    composition = estimate_ancestry(genotypes)
    admixture = detect_admixture(genotypes)
    
    return {
        "disclaimer": ANCESTRY_DISCLAIMER,
        "composition": composition,
        "admixture": admixture,
        "summary": _generate_ancestry_text_summary(composition, admixture),
        "methodology": {
            "approach": "Likelihood-based estimation using ancestry informative markers",
            "resolution": "Continental level only",
            "limitations": [
                "Sub-regional ancestry (countries, ethnic groups) cannot be determined",
                "All estimates include wide confidence intervals",
                "Reference panels may not represent all populations equally"
            ],
            "pmid": ["26432245", "12493913", "18292342"]
        }
    }


def _generate_ancestry_text_summary(composition: Dict, admixture: Dict) -> str:
    """Generate text summary for reports."""
    lines = []
    
    lines.append("⚠️ CONTINENTAL-LEVEL ESTIMATES ONLY")
    lines.append("Sub-regional ancestry (e.g., specific countries) cannot be reliably determined.")
    lines.append("")
    
    if composition.get("status") != "success":
        return "Insufficient ancestry informative markers for reliable estimation."
    
    primary = composition['primary_ancestry']
    pct = composition['primary_percentage']
    range_str = composition.get('primary_range', 'unknown')
    
    lines.append(f"Primary ancestry: {primary} — {pct}% (range: {range_str})")
    lines.append("")
    
    lines.append("Full breakdown with confidence intervals:")
    for pop, display in composition.get("ancestry_proportions", {}).items():
        lines.append(f"  {pop}: {display}")
    
    lines.append("")
    
    if admixture.get("is_admixed"):
        lines.append("Admixture: " + admixture.get("interpretation", "Detected"))
    else:
        lines.append("Admixture: No significant multi-continental admixture detected")
    
    lines.append("")
    lines.append(f"Confidence: {composition.get('confidence', 'unknown')}")
    lines.append(f"Note: {composition.get('confidence_note', '')}")
    lines.append(f"Markers analyzed: {composition.get('markers_used', 0)}")
    
    return "\n".join(lines)
