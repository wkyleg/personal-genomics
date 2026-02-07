"""
Ancient Ancestral Signals Module
Detection of ancestral population signatures in modern genomes

This module detects signals from ancient ancestral populations:
- Western Hunter-Gatherers (WHG) - Mesolithic Europeans (~15,000-8,000 BP)
- Early European Farmers (EEF) - Anatolian Neolithic (~10,000-5,000 BP)
- Steppe Pastoralists - Yamnaya/Pontic-Caspian (~5,000-4,000 BP)
- Archaic Introgression - Neanderthal (~50,000-40,000 BP) and Denisovan

╔══════════════════════════════════════════════════════════════════════════════╗
║  WHY ANCIENT SIGNALS INSTEAD OF MODERN ETHNICITY?                            ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Modern ethnicity percentages (e.g., "32% Irish") are scientifically         ║
║  problematic because:                                                         ║
║                                                                              ║
║  • Consumer arrays can't distinguish closely related populations              ║
║  • Modern national boundaries don't reflect genetic history                   ║
║  • Ancient migrations mixed populations extensively                           ║
║                                                                              ║
║  Ancient ancestral signals are:                                              ║
║  • Based on well-characterized ancient DNA sequences                          ║
║  • Scientifically validated with clear trait associations                     ║
║  • More informative about actual genetic heritage                             ║
╚══════════════════════════════════════════════════════════════════════════════╝

Key References:
- Lazaridis et al. 2014. Ancient human genomes suggest three ancestral 
  populations for present-day Europeans. PMID: 25230663
- Haak et al. 2015. Massive migration from the steppe was a source for 
  Indo-European languages in Europe. PMID: 25731166
- Mathieson et al. 2015. Genome-wide patterns of selection in 230 ancient 
  Eurasians. PMID: 26595274
- Prüfer et al. 2014. The complete genome sequence of a Neanderthal. 
  PMID: 24352235
- Vernot & Akey 2014. Resurrecting surviving Neandertal lineages from 
  modern human genomes. PMID: 24476815
"""

from typing import Dict, List, Optional, Any, Tuple
import math

# =============================================================================
# ANCIENT POPULATION DESCRIPTIONS
# =============================================================================

ANCIENT_POPULATIONS = {
    "WHG": {
        "name": "Western Hunter-Gatherers (Mesolithic)",
        "abbreviation": "WHG",
        "time_period": "~15,000-8,000 years ago",
        "description": "Indigenous European hunter-gatherer population before the arrival of farming. "
                      "Survived in Mesolithic Europe, later mixed with incoming farmers.",
        "geographic_origin": "Western and Central Europe (Loschbour, La Braña, Cheddar Man)",
        "traits_contributed": [
            "Blue eyes (likely origin of European blue eyes)",
            "Dark skin (contrary to modern stereotypes)",
            "Lactose intolerance (original European state)",
            "Immune adaptations to European pathogens"
        ],
        "modern_distribution": "Present in all Europeans at varying levels (~10-30%). "
                              "Highest in Baltic and Scandinavian populations.",
        "pmid": ["25230663", "24463515", "26595274"]
    },
    "EEF": {
        "name": "Early European Farmers (Neolithic Anatolian)",
        "abbreviation": "EEF",
        "time_period": "~10,000-5,000 years ago",
        "description": "Anatolian farmers who migrated to Europe bringing agriculture. "
                      "Largely replaced WHG populations but mixed with them over time.",
        "geographic_origin": "Anatolia (modern Turkey), spread through Mediterranean and Central Europe",
        "traits_contributed": [
            "Lighter skin (adaptation to lower vitamin D from grain diet)",
            "Brown eyes (ancestral state)",
            "Shorter stature than WHG",
            "Metabolic adaptations to agricultural diet",
            "Some lactose tolerance (not as strong as later Steppe)"
        ],
        "modern_distribution": "Major component of Southern European ancestry (~60-80%). "
                              "Present throughout Europe, decreasing northward.",
        "pmid": ["25230663", "25731166", "26595274"]
    },
    "STEPPE": {
        "name": "Steppe Pastoralists (Yamnaya/Pontic-Caspian)",
        "abbreviation": "Steppe",
        "time_period": "~5,000-4,000 years ago",
        "description": "Bronze Age pastoralists from the Pontic-Caspian steppe who migrated "
                      "into Europe. Brought Indo-European languages, horses, and the wheel.",
        "geographic_origin": "Pontic-Caspian Steppe (Ukraine/Southern Russia)",
        "traits_contributed": [
            "Strong lactose tolerance (LCT persistence)",
            "Taller stature",
            "Lighter skin and hair pigmentation",
            "Indo-European language family",
            "Some immunity genes"
        ],
        "modern_distribution": "Major component of Northern/Eastern European ancestry (~30-50%). "
                              "Highest in Northern Europe, also present in South Asia (Indo-Aryan).",
        "pmid": ["25731166", "25230663", "26062507"]
    },
    "NEANDERTHAL": {
        "name": "Neanderthal Introgression",
        "abbreviation": "Neanderthal",
        "time_period": "~50,000-40,000 years ago (introgression events)",
        "description": "Archaic human species that interbred with anatomically modern humans "
                      "as they left Africa. All non-African humans carry ~1-4% Neanderthal DNA.",
        "geographic_origin": "Europe and Western Asia",
        "traits_contributed": [
            "Immune system genes (HLA variants)",
            "Keratin/hair/skin genes",
            "Some pain sensitivity variants",
            "Lipid metabolism genes",
            "Some linked to depression, blood clotting"
        ],
        "modern_distribution": "~1-4% in all non-Africans. Europeans and East Asians carry "
                              "different Neanderthal segments. Sub-Saharan Africans have minimal.",
        "pmid": ["24352235", "24476815", "26194313", "29058716"]
    },
    "DENISOVAN": {
        "name": "Denisovan Introgression",
        "abbreviation": "Denisovan",
        "time_period": "~50,000-30,000 years ago (introgression events)",
        "description": "Archaic human species known primarily from DNA. Interbred with ancestors "
                      "of Melanesians, Aboriginal Australians, and some East/South Asians.",
        "geographic_origin": "Siberia/Central Asia (Denisova Cave), likely broader Asian range",
        "traits_contributed": [
            "High-altitude adaptation (EPAS1 in Tibetans)",
            "Immune genes",
            "Fat metabolism",
            "Possible cold adaptation"
        ],
        "modern_distribution": "Highest in Melanesians/Papuans (~3-6%). Lower in East Asians, "
                              "South Asians. Minimal in Europeans and Africans.",
        "pmid": ["24352235", "25043035", "29058716"]
    }
}

# =============================================================================
# ANCIENT ANCESTRY MARKERS
# =============================================================================

# Markers with strong ancient population associations
# Based on ancient DNA studies and selection scans

ANCIENT_ANCESTRY_MARKERS = {
    # =========================================================================
    # WESTERN HUNTER-GATHERER (WHG) SIGNALS
    # =========================================================================
    
    # Blue eyes - WHG signature (later spread by Steppe)
    "rs12913832": {
        "gene": "HERC2/OCA2",
        "trait": "Blue eyes",
        "ancestral_population": "WHG",
        "derived_allele": "G",
        "ancestral_allele": "A",
        "population_frequencies": {
            "WHG": 0.90,  # Ancient samples show high frequency
            "EEF": 0.05,  # Rare in early farmers
            "Steppe": 0.50,  # Intermediate
        },
        "signal_interpretation": {
            "GG": "Strong WHG signal (blue eyes likely originated in Mesolithic Europe)",
            "AG": "Moderate WHG signal",
            "AA": "Ancestral state (brown eyes, minimal WHG signal at this locus)"
        },
        "note": "Blue eyes were present in WHG like Loschbour and La Braña individuals",
        "pmid": ["24463515", "25230663", "26595274"]
    },
    
    # SLC45A2 - Skin pigmentation (WHG had dark skin!)
    "rs16891982": {
        "gene": "SLC45A2",
        "trait": "Skin pigmentation",
        "ancestral_population": "EEF/Steppe",  # Light skin came LATER
        "derived_allele": "G",
        "ancestral_allele": "C",
        "population_frequencies": {
            "WHG": 0.05,  # WHG had dark skin
            "EEF": 0.80,  # Farmers brought lighter skin
            "Steppe": 0.95,  # Strongly selected
        },
        "signal_interpretation": {
            "GG": "Derived state (lighter skin - EEF/Steppe signal)",
            "CG": "Heterozygous",
            "CC": "Ancestral state (darker skin - possible WHG retention)"
        },
        "note": "Contrary to stereotypes, European hunter-gatherers had dark skin",
        "pmid": ["25230663", "26595274", "24463515"]
    },
    
    # =========================================================================
    # EARLY EUROPEAN FARMER (EEF) SIGNALS
    # =========================================================================
    
    # SLC24A5 - Major skin lightening (selected in farmers)
    "rs1426654": {
        "gene": "SLC24A5",
        "trait": "Skin pigmentation",
        "ancestral_population": "EEF",
        "derived_allele": "A",
        "ancestral_allele": "G",
        "population_frequencies": {
            "WHG": 0.10,  # Low in hunter-gatherers
            "EEF": 0.95,  # High in Neolithic farmers
            "Steppe": 0.95,  # Also high
        },
        "signal_interpretation": {
            "AA": "Strong EEF/post-Neolithic signal (nearly fixed in Europeans)",
            "AG": "Heterozygous (rare in Europeans)",
            "GG": "Ancestral African/WHG-like state (very rare in Europeans)"
        },
        "note": "Selected for vitamin D synthesis on grain-based diet with less fish",
        "pmid": ["16357253", "26595274", "25230663"]
    },
    
    # NAT2 - Metabolic adaptation to plant-based diet
    "rs1801280": {
        "gene": "NAT2",
        "trait": "Acetylator status (drug/food metabolism)",
        "ancestral_population": "EEF",
        "derived_allele": "A",
        "ancestral_allele": "G",
        "population_frequencies": {
            "WHG": 0.20,
            "EEF": 0.60,
            "Steppe": 0.40,
        },
        "signal_interpretation": {
            "AA": "Slow acetylator (EEF dietary adaptation signal)",
            "AG": "Intermediate",
            "GG": "Fast acetylator"
        },
        "note": "Slow acetylation may have been advantageous for plant toxin metabolism",
        "pmid": ["26595274", "25731166"]
    },
    
    # =========================================================================
    # STEPPE PASTORALIST SIGNALS
    # =========================================================================
    
    # LCT - Lactase persistence (STRONG Steppe signature)
    "rs4988235": {
        "gene": "LCT",
        "trait": "Lactase persistence (adult milk digestion)",
        "ancestral_population": "Steppe",
        "derived_allele": "A",
        "ancestral_allele": "G",
        "population_frequencies": {
            "WHG": 0.00,  # Absent in hunter-gatherers
            "EEF": 0.05,  # Very rare in early farmers
            "Steppe": 0.30,  # Present, then rapidly selected
        },
        "signal_interpretation": {
            "AA": "Strong lactase persistence (Steppe pastoralist signature)",
            "AG": "Likely lactase persistent (Steppe signal)",
            "GG": "Lactase non-persistent (ancestral pre-Steppe state)"
        },
        "note": "Rose from ~10% to ~80% in Europe in just 4,000 years - strongest "
               "known selection in recent human evolution. Steppe pastoralists "
               "brought cattle-herding culture.",
        "pmid": ["25731166", "26595274", "26062507"]
    },
    
    # FADS genes - Fatty acid metabolism (Steppe diet adaptation)
    "rs174546": {
        "gene": "FADS1",
        "trait": "Fatty acid desaturation",
        "ancestral_population": "Steppe",
        "derived_allele": "T",
        "ancestral_allele": "C",
        "population_frequencies": {
            "WHG": 0.30,
            "EEF": 0.40,
            "Steppe": 0.70,
        },
        "signal_interpretation": {
            "TT": "Efficient plant-based omega-3 conversion (selected post-Steppe)",
            "CT": "Intermediate",
            "CC": "Relies more on dietary omega-3 (fish-eating adaptation)"
        },
        "note": "Selection for converting plant omega-3 to brain-essential DHA",
        "pmid": ["26595274", "27182965"]
    },
    
    # Height-associated variant (Steppe were taller)
    "rs1042725": {
        "gene": "HMGA2",
        "trait": "Height",
        "ancestral_population": "Steppe",
        "derived_allele": "C",
        "ancestral_allele": "T",
        "population_frequencies": {
            "WHG": 0.40,
            "EEF": 0.30,  # Farmers were shorter
            "Steppe": 0.60,  # Steppe were tall
        },
        "signal_interpretation": {
            "CC": "Height-increasing allele (Steppe signal)",
            "CT": "Intermediate",
            "TT": "Height-decreasing allele"
        },
        "note": "Steppe pastoralists were significantly taller than Neolithic farmers",
        "pmid": ["26595274", "25731166"]
    },
    
    # =========================================================================
    # NEANDERTHAL INTROGRESSION MARKERS
    # =========================================================================
    
    # BNC2 - Skin/freckling (Neanderthal origin)
    "rs10756819": {
        "gene": "BNC2",
        "trait": "Skin pigmentation, freckling",
        "ancestral_population": "Neanderthal",
        "derived_allele": "G",
        "ancestral_allele": "A",
        "population_frequencies": {
            "Neanderthal": 0.90,
            "Modern_EUR": 0.70,
            "Modern_AFR": 0.05,
        },
        "signal_interpretation": {
            "GG": "Neanderthal-introgressed allele (associated with lighter skin, freckling)",
            "AG": "Heterozygous for Neanderthal variant",
            "AA": "Modern human ancestral state"
        },
        "note": "This region shows strong evidence of adaptive introgression from Neanderthals",
        "pmid": ["24476815", "26194313", "28133863"]
    },
    
    # OAS1 - Immune function (Neanderthal-derived)
    "rs10774671": {
        "gene": "OAS1",
        "trait": "Antiviral immune response",
        "ancestral_population": "Neanderthal",
        "derived_allele": "G",
        "ancestral_allele": "A",
        "population_frequencies": {
            "Neanderthal": 1.00,
            "Modern_EUR": 0.35,
            "Modern_EAS": 0.40,
            "Modern_AFR": 0.05,
        },
        "signal_interpretation": {
            "GG": "Neanderthal-introgressed immune variant",
            "AG": "Heterozygous",
            "AA": "Modern human ancestral state"
        },
        "note": "Neanderthal OAS variants may provide enhanced antiviral defense",
        "pmid": ["24476815", "26194313", "27654910"]
    },
    
    # TLR genes - Pathogen response (Neanderthal origin)
    "rs5743810": {
        "gene": "TLR6",
        "trait": "Pathogen recognition, immune response",
        "ancestral_population": "Neanderthal",
        "derived_allele": "G",
        "ancestral_allele": "A",
        "population_frequencies": {
            "Neanderthal": 0.95,
            "Modern_EUR": 0.45,
            "Modern_EAS": 0.30,
            "Modern_AFR": 0.02,
        },
        "signal_interpretation": {
            "GG": "Neanderthal-introgressed immune variant",
            "AG": "Heterozygous for Neanderthal variant",
            "AA": "Modern human ancestral state"
        },
        "note": "Toll-like receptor variants from Neanderthals help fight pathogens",
        "pmid": ["26194313", "29058716"]
    },
    
    # SLC16A11 - Type 2 diabetes risk (Neanderthal)
    "rs13342232": {
        "gene": "SLC16A11",
        "trait": "Type 2 diabetes risk, lipid metabolism",
        "ancestral_population": "Neanderthal",
        "derived_allele": "T",
        "ancestral_allele": "C",
        "population_frequencies": {
            "Neanderthal": 1.00,
            "Modern_AMR": 0.50,  # High in Native Americans/Mexicans
            "Modern_EAS": 0.10,
            "Modern_EUR": 0.02,
            "Modern_AFR": 0.00,
        },
        "signal_interpretation": {
            "TT": "Neanderthal-derived diabetes risk variant",
            "CT": "Heterozygous carrier",
            "CC": "Modern human ancestral (lower risk)"
        },
        "note": "Neanderthal variant that increases T2D risk, common in Native Americans",
        "pmid": ["24476815", "25282103"]
    },
    
    # =========================================================================
    # DENISOVAN MARKERS
    # =========================================================================
    
    # EPAS1 - High altitude adaptation (Denisovan in Tibetans)
    "rs115321619": {
        "gene": "EPAS1",
        "trait": "High altitude adaptation",
        "ancestral_population": "Denisovan",
        "derived_allele": "G",
        "ancestral_allele": "A",
        "population_frequencies": {
            "Denisovan": 1.00,
            "Tibetan": 0.87,  # Strong selection at altitude
            "Han_Chinese": 0.09,
            "Modern_EUR": 0.00,
        },
        "signal_interpretation": {
            "GG": "Denisovan-introgressed high-altitude allele (Tibetan signature)",
            "AG": "Heterozygous (rare outside Tibet)",
            "AA": "Normal low-altitude state"
        },
        "note": "Denisovan variant enables Tibetans to thrive at high altitude. "
               "Classic example of adaptive introgression.",
        "pmid": ["25043035", "24976165"]
    },
    
    # TBX15/WARS2 - Body fat distribution (possible Denisovan)
    "rs2298080": {
        "gene": "TBX15",
        "trait": "Body fat distribution, ear morphology",
        "ancestral_population": "Denisovan",
        "derived_allele": "C",
        "ancestral_allele": "T",
        "population_frequencies": {
            "Denisovan": 0.95,
            "Modern_EAS": 0.65,
            "Inuit": 0.80,
            "Modern_EUR": 0.20,
        },
        "signal_interpretation": {
            "CC": "Possible Denisovan-introgressed body fat/cold adaptation",
            "CT": "Heterozygous",
            "TT": "Modern human ancestral state"
        },
        "note": "May be related to cold adaptation, high in Arctic populations",
        "pmid": ["25043035", "29058716"]
    },
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def detect_ancient_signals(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Detect ancient ancestral population signals in genotype data.
    
    Returns detected signals (not percentages) for each ancient population,
    along with which specific markers were found.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Dict with signals detected for each ancient population
    """
    results = {
        "methodology": {
            "description": "Detection of ancient ancestral population signals using "
                          "ancestry-informative markers from ancient DNA studies",
            "approach": "Signal detection, not percentage estimation",
            "note": "Signals indicate presence of ancestry; strength reflects marker count, "
                   "not precise ancestry proportions",
            "pmid": ["25230663", "25731166", "24352235"]
        },
        "populations": {},
        "markers_found": 0,
        "markers_checked": len(ANCIENT_ANCESTRY_MARKERS)
    }
    
    # Initialize population signal tracking
    population_signals = {
        "WHG": {"markers_found": [], "signal_strength": "none"},
        "EEF": {"markers_found": [], "signal_strength": "none"},
        "STEPPE": {"markers_found": [], "signal_strength": "none"},
        "NEANDERTHAL": {"markers_found": [], "signal_strength": "none"},
        "DENISOVAN": {"markers_found": [], "signal_strength": "none"},
    }
    
    # Check each marker
    for rsid, info in ANCIENT_ANCESTRY_MARKERS.items():
        geno = genotypes.get(rsid)
        if not geno:
            continue
            
        results["markers_found"] += 1
        
        # Determine which population this marker signals
        anc_pop = info.get("ancestral_population", "").upper()
        if anc_pop == "WHG":
            pop_key = "WHG"
        elif anc_pop in ["EEF", "NEOLITHIC"]:
            pop_key = "EEF"
        elif anc_pop in ["STEPPE", "YAMNAYA"]:
            pop_key = "STEPPE"
        elif anc_pop == "NEANDERTHAL":
            pop_key = "NEANDERTHAL"
        elif anc_pop == "DENISOVAN":
            pop_key = "DENISOVAN"
        else:
            continue
        
        # Check if derived allele is present
        derived = info.get("derived_allele", "")
        derived_count = geno.upper().count(derived.upper()) if derived else 0
        
        # Get interpretation
        interp = info.get("signal_interpretation", {}).get(geno, "")
        
        if derived_count > 0:
            population_signals[pop_key]["markers_found"].append({
                "rsid": rsid,
                "gene": info.get("gene", ""),
                "trait": info.get("trait", ""),
                "genotype": geno,
                "derived_copies": derived_count,
                "interpretation": interp,
                "pmid": info.get("pmid", [])
            })
    
    # Calculate signal strength for each population
    for pop_key, data in population_signals.items():
        marker_count = len(data["markers_found"])
        
        if marker_count == 0:
            data["signal_strength"] = "not detected"
            data["signal_description"] = "No markers found for this ancestral population"
        elif marker_count == 1:
            data["signal_strength"] = "weak"
            data["signal_description"] = "Single marker detected - signal present but limited"
        elif marker_count == 2:
            data["signal_strength"] = "moderate"
            data["signal_description"] = "Multiple markers detected - clear signal present"
        else:
            data["signal_strength"] = "strong"
            data["signal_description"] = "Multiple markers detected - strong ancestral signal"
        
        # Add population info
        pop_info = ANCIENT_POPULATIONS.get(pop_key, {})
        results["populations"][pop_key] = {
            "name": pop_info.get("name", pop_key),
            "time_period": pop_info.get("time_period", "Unknown"),
            "description": pop_info.get("description", ""),
            "traits_contributed": pop_info.get("traits_contributed", []),
            "signal_strength": data["signal_strength"],
            "signal_description": data["signal_description"],
            "markers_detected": data["markers_found"],
            "marker_count": marker_count,
            "pmid": pop_info.get("pmid", [])
        }
    
    return results


def get_ancestry_summary(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Complete ancient ancestry analysis.
    
    Returns signal detection results for ancient populations,
    NOT modern ethnicity percentages.
    """
    signals = detect_ancient_signals(genotypes)
    
    return {
        "analysis_type": "ancient_ancestral_signals",
        "disclaimer": (
            "This analysis detects SIGNALS from ancient ancestral populations, "
            "not modern ethnicity percentages. Modern ethnicity estimates are "
            "scientifically unreliable from consumer DNA arrays. Ancient signals "
            "are based on well-characterized ancient DNA studies."
        ),
        "signals": signals,
        "summary": _generate_signals_summary(signals),
        "educational_note": (
            "All non-African humans are ~98% 'anatomically modern human' with "
            "~1-4% archaic (Neanderthal/Denisovan) ancestry. Europeans are a "
            "mixture of three main ancient populations: Mesolithic hunter-gatherers, "
            "Neolithic farmers from Anatolia, and Bronze Age Steppe pastoralists."
        )
    }


def _generate_signals_summary(signals: Dict) -> str:
    """Generate human-readable summary of ancient signals."""
    lines = []
    
    lines.append("ANCIENT ANCESTRAL SIGNALS DETECTED")
    lines.append("=" * 40)
    lines.append("")
    lines.append("Note: These are signal detections, not precise percentages.")
    lines.append("Signal strength reflects marker availability, not ancestry proportion.")
    lines.append("")
    
    for pop_key in ["WHG", "EEF", "STEPPE", "NEANDERTHAL", "DENISOVAN"]:
        pop = signals.get("populations", {}).get(pop_key, {})
        name = pop.get("name", pop_key)
        strength = pop.get("signal_strength", "unknown")
        period = pop.get("time_period", "")
        count = pop.get("marker_count", 0)
        
        # Signal indicator
        if strength == "strong":
            indicator = "●●●"
        elif strength == "moderate":
            indicator = "●●○"
        elif strength == "weak":
            indicator = "●○○"
        else:
            indicator = "○○○"
        
        lines.append(f"{name}")
        lines.append(f"  Time: {period}")
        lines.append(f"  Signal: {indicator} {strength.upper()} ({count} markers)")
        
        if count > 0:
            traits = pop.get("traits_contributed", [])[:3]
            if traits:
                lines.append(f"  Traits: {', '.join(traits)}")
        
        lines.append("")
    
    return "\n".join(lines)


# =============================================================================
# LEGACY COMPATIBILITY
# =============================================================================

# Keep old function names for compatibility, but they now use ancient signals

def estimate_ancestry(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Legacy compatibility function.
    Now returns ancient signals instead of modern ethnicity estimates.
    """
    return get_ancestry_summary(genotypes)


def detect_admixture(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Legacy compatibility function.
    Admixture in ancient context means mixing of ancient populations.
    """
    signals = detect_ancient_signals(genotypes)
    
    # Count populations with signals
    detected_pops = [
        p for p, data in signals.get("populations", {}).items()
        if data.get("signal_strength") not in ["not detected", "none"]
    ]
    
    return {
        "analysis_type": "ancient_population_mixture",
        "note": "All modern Europeans are mixtures of ancient populations",
        "populations_detected": detected_pops,
        "interpretation": (
            f"Signals from {len(detected_pops)} ancient populations detected. "
            "This is expected - modern humans outside Africa are mixtures of "
            "ancient populations including archaic humans (Neanderthal/Denisovan)."
        )
    }


# Wilson CI function kept for PRS module
def calculate_wilson_confidence_interval(
    successes: int,
    trials: int,
    confidence: float = 0.95
) -> Tuple[float, float]:
    """
    Calculate Wilson score confidence interval for a proportion.
    Kept for use by other modules (PRS).
    """
    if trials == 0:
        return (0.0, 1.0)
    
    z = 1.96 if confidence == 0.95 else 1.645 if confidence == 0.90 else 2.576
    
    p_hat = successes / trials
    
    denominator = 1 + z**2 / trials
    center = (p_hat + z**2 / (2 * trials)) / denominator
    
    margin = z * math.sqrt((p_hat * (1 - p_hat) + z**2 / (4 * trials)) / trials) / denominator
    
    lower = max(0, center - margin)
    upper = min(1, center + margin)
    
    return (lower, upper)
