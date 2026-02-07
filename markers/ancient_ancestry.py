"""
Ancient DNA Analysis Module

Detects genetic signals from ancient human populations:
- Western Hunter-Gatherers (WHG) - Mesolithic Europe
- Eastern Hunter-Gatherers (EHG) - Mesolithic Russia
- Anatolian Neolithic Farmers (ANF) - Early farmers
- Yamnaya/Steppe Pastoralists - Bronze Age herders
- Neanderthal introgression markers
- Denisovan introgression markers (limited detection on arrays)

This module provides EDUCATIONAL context about the deep history of human migrations,
not precise ancestry percentages.

Key Sources:
- Lazaridis et al. 2014 - Ancient human genomes suggest three ancestral populations (PMID: 25230663)
- Haak et al. 2015 - Massive migration from the steppe (PMID: 25731166)
- Mathieson et al. 2015 - Genome-wide patterns of selection in Europe (PMID: 26595274)
- Prüfer et al. 2014 - The complete genome sequence of a Neanderthal (PMID: 24352235)
"""

from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
import json

# =============================================================================
# LOAD REFERENCE DATA
# =============================================================================

def load_ancient_dna_markers() -> Dict[str, Any]:
    """Load the ancient DNA markers reference data."""
    ref_path = Path(__file__).parent.parent / "references" / "ancient_dna_markers.json"
    
    if not ref_path.exists():
        return {"_metadata": {"ancient_populations": {}}}
    
    with open(ref_path, 'r') as f:
        return json.load(f)


# Cache the reference data
_ANCIENT_DATA = None

def get_ancient_data() -> Dict[str, Any]:
    """Get cached ancient DNA reference data."""
    global _ANCIENT_DATA
    if _ANCIENT_DATA is None:
        _ANCIENT_DATA = load_ancient_dna_markers()
    return _ANCIENT_DATA


# =============================================================================
# ANCIENT POPULATION METADATA
# =============================================================================

ANCIENT_POPULATIONS = {
    "WHG": {
        "name": "Western Hunter-Gatherers",
        "emoji": "🏹",
        "period": "~15,000-5,000 BCE",
        "region": "Western & Central Europe",
        "description": "The original inhabitants of post-Ice Age Europe. They were foragers who lived in small groups across the continent for thousands of years before farmers arrived.",
        "physical": "Dark skin, blue eyes (a surprising combination not seen today). This has been confirmed by ancient DNA from sites like Loschbour (Luxembourg) and La Braña (Spain).",
        "key_sites": ["Loschbour (Luxembourg)", "La Braña (Spain)", "Cheddar Man (UK)"],
        "modern_contribution": "10-20% of modern European ancestry",
        "pmid": ["25230663", "25731166"]
    },
    "EHG": {
        "name": "Eastern Hunter-Gatherers",
        "emoji": "🦌",
        "period": "~8,000-5,000 BCE",
        "region": "Eastern Europe & Russia",
        "description": "Hunter-gatherers from the forest-steppe zone extending from the Baltic to Siberia. They were ancestral to later Steppe populations.",
        "physical": "Mixed pigmentation - some had lighter skin than WHG. Closer genetic ties to Ancient North Eurasians.",
        "key_sites": ["Karelia (Russia)", "Samara (Russia)"],
        "modern_contribution": "Part of Steppe ancestry in Europeans",
        "pmid": ["25230663", "25731166"]
    },
    "ANF": {
        "name": "Anatolian Neolithic Farmers",
        "emoji": "🌾",
        "period": "~8,000-4,000 BCE",
        "region": "Anatolia (Turkey) → Europe",
        "description": "The first farmers who migrated from the Near East into Europe starting ~9,000 years ago. They brought agriculture, domesticated animals, and new ways of life.",
        "physical": "Light skin (had the modern European SLC24A5 variant), brown eyes. Shorter stature than hunter-gatherers.",
        "key_sites": ["Barcın (Turkey)", "Stuttgart (Germany)", "Otzi the Iceman"],
        "modern_contribution": "30-50% of modern European ancestry, highest in Mediterranean",
        "pmid": ["25230663", "26595274"]
    },
    "Yamnaya": {
        "name": "Yamnaya/Steppe Pastoralists",
        "emoji": "🐴",
        "period": "~3,300-2,600 BCE",
        "region": "Pontic-Caspian Steppe",
        "description": "Bronze Age herders who domesticated the horse and invented wheeled vehicles. They migrated massively into Europe ~5,000 years ago, transforming its genetic landscape. Likely spread Indo-European languages.",
        "physical": "Light skin, variable eye color (mix of WHG and ANF variants), tall stature.",
        "key_sites": ["Samara (Russia)", "Yamnaya kurgans (Ukraine/Russia)"],
        "modern_contribution": "30-50% of modern Northern European ancestry",
        "pmid": ["25731166", "26595274"]
    },
    "Neanderthal": {
        "name": "Neanderthal",
        "emoji": "🦴",
        "period": "~400,000-40,000 years ago",
        "region": "Europe & Western Asia",
        "description": "Our archaic human relatives who lived in Europe for hundreds of thousands of years. Modern humans interbred with them after leaving Africa, leaving ~1-4% Neanderthal DNA in all non-Africans.",
        "physical": "Robust build, large nose (possibly cold adaptation), larger brain than modern humans.",
        "key_sites": ["Vindija Cave (Croatia)", "Altai Mountains (Russia)", "El Sidrón (Spain)"],
        "modern_contribution": "~1-4% in non-African populations (0% in sub-Saharan Africans)",
        "pmid": ["24352235", "23128226"]
    },
    "Denisovan": {
        "name": "Denisovan",
        "emoji": "🏔️",
        "period": "~100,000-30,000 years ago",
        "region": "Central & East Asia",
        "description": "A mysterious archaic human group known mainly from DNA. Related to Neanderthals but distinct. Modern humans interbred with them too, especially in East Asia and Oceania.",
        "physical": "Largely unknown - only a few bone fragments exist. May have been adapted to high altitudes.",
        "key_sites": ["Denisova Cave (Siberia)"],
        "modern_contribution": "~3-5% in Melanesians, <1% in East Asians",
        "pmid": ["21179161", "24352235"]
    }
}


def normalize_genotype(genotype: str) -> str:
    """Normalize genotype to alphabetical order."""
    if len(genotype) == 2:
        return "".join(sorted(genotype.upper()))
    return genotype.upper()


# =============================================================================
# ANCIENT SIGNAL DETECTION
# =============================================================================

def detect_ancient_signals(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Detect markers associated with ancient populations.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Dict with detected signals for each ancient population
    """
    ancient_data = get_ancient_data()
    
    results = {
        "WHG": {"markers_detected": [], "total_markers": 0},
        "EHG": {"markers_detected": [], "total_markers": 0},
        "ANF": {"markers_detected": [], "total_markers": 0},
        "Yamnaya": {"markers_detected": [], "total_markers": 0},
        "Neanderthal": {"markers_detected": [], "total_markers": 0},
        "Denisovan": {"markers_detected": [], "total_markers": 0},
    }
    
    for rsid, marker_data in ancient_data.items():
        if rsid.startswith("_"):
            continue
        
        user_geno = genotypes.get(rsid)
        if not user_geno:
            continue
        
        user_geno_norm = normalize_genotype(user_geno)
        ancient_context = marker_data.get("ancient_context", {})
        introgression = marker_data.get("introgression", {})
        
        # Check each ancient population
        for pop_key, pop_data in ancient_context.items():
            # Map alternative keys to our standard population names
            pop_mapped = {
                "mesolithic_whg": "WHG",
                "neolithic_farmer": "ANF",
                "yamnaya_steppe": "Yamnaya",
                "AncientEastAsian": "EHG",  # Close enough for display
                "AncientAfrican": None,  # Skip - not in our ancient population list
            }.get(pop_key, pop_key)
            
            if pop_mapped not in results:
                continue
            
            results[pop_mapped]["total_markers"] += 1
            
            genotypes_found = pop_data.get("genotypes_found", [])
            freq = pop_data.get("frequency", "unknown")
            note = pop_data.get("note", "")
            
            # Check if user's genotype matches ancient pattern
            for ancient_geno in genotypes_found:
                if normalize_genotype(ancient_geno) == user_geno_norm:
                    results[pop_mapped]["markers_detected"].append({
                        "rsid": rsid,
                        "gene": marker_data.get("gene", ""),
                        "trait": marker_data.get("trait", ""),
                        "your_genotype": user_geno,
                        "ancient_frequency": freq,
                        "note": note,
                        "history": marker_data.get("history", ""),
                        "pmid": marker_data.get("pmid", [])
                    })
                    break
        
        # Check for Neanderthal/Denisovan introgression
        if introgression:
            source = introgression.get("source", "")
            if source in results:
                results[source]["total_markers"] += 1
                
                derived = marker_data.get("derived_allele", "")
                if derived and derived.upper() in user_geno.upper():
                    results[source]["markers_detected"].append({
                        "rsid": rsid,
                        "gene": marker_data.get("gene", ""),
                        "trait": marker_data.get("trait", ""),
                        "your_genotype": user_geno,
                        "introgression_source": source,
                        "modern_frequency_EUR": introgression.get("modern_frequency_EUR", "unknown"),
                        "effect": introgression.get("effect", ""),
                        "history": marker_data.get("history", ""),
                        "pmid": marker_data.get("pmid", [])
                    })
    
    # Calculate detection summaries
    for pop in results:
        markers = results[pop]["markers_detected"]
        total = results[pop]["total_markers"]
        results[pop]["detection_count"] = len(markers)
        results[pop]["detection_rate"] = len(markers) / total if total > 0 else 0
        results[pop]["metadata"] = ANCIENT_POPULATIONS.get(pop, {})
    
    return results


def generate_ancient_dna_report(genotypes: Dict[str, str]) -> str:
    """
    Generate a narrative report about ancient ancestral signals.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Formatted text report
    """
    signals = detect_ancient_signals(genotypes)
    
    lines = []
    lines.append("=" * 70)
    lines.append("🏛️ ANCIENT ANCESTRAL SIGNALS")
    lines.append("=" * 70)
    lines.append("")
    lines.append("This analysis examines your DNA for markers associated with ancient")
    lines.append("human populations. These are EDUCATIONAL insights about deep history,")
    lines.append("not precise ancestry percentages.")
    lines.append("")
    lines.append("⚠️  Consumer DNA arrays can only detect a handful of ancient markers.")
    lines.append("Full ancient DNA analysis requires whole genome sequencing.")
    lines.append("")
    
    # Order populations chronologically
    pop_order = ["WHG", "EHG", "ANF", "Yamnaya", "Neanderthal", "Denisovan"]
    
    for pop_key in pop_order:
        pop_data = signals.get(pop_key, {})
        metadata = ANCIENT_POPULATIONS.get(pop_key, {})
        
        emoji = metadata.get("emoji", "")
        name = metadata.get("name", pop_key)
        period = metadata.get("period", "")
        
        lines.append("-" * 70)
        lines.append(f"{emoji} {name.upper()}")
        lines.append(f"   {period}")
        lines.append("-" * 70)
        lines.append("")
        
        # Description
        description = metadata.get("description", "")
        if description:
            # Word wrap
            words = description.split()
            line = "   "
            for word in words:
                if len(line) + len(word) > 65:
                    lines.append(line)
                    line = "   "
                line += word + " "
            if line.strip():
                lines.append(line)
            lines.append("")
        
        # Physical traits
        physical = metadata.get("physical", "")
        if physical:
            lines.append(f"   Physical: {physical}")
            lines.append("")
        
        # Detected markers
        markers = pop_data.get("markers_detected", [])
        
        if markers:
            lines.append("   DETECTED MARKERS:")
            for marker in markers:
                rsid = marker.get("rsid", "")
                gene = marker.get("gene", "")
                geno = marker.get("your_genotype", "")
                trait = marker.get("trait", "")
                note = marker.get("note", "")
                
                gene_str = f" ({gene})" if gene else ""
                lines.append(f"   ✓ {rsid}{gene_str}: {geno}")
                if trait:
                    lines.append(f"      Trait: {trait}")
                if note:
                    # Word wrap note
                    note_words = note.split()
                    note_line = "      → "
                    for word in note_words:
                        if len(note_line) + len(word) > 60:
                            lines.append(note_line)
                            note_line = "        "
                        note_line += word + " "
                    if note_line.strip():
                        lines.append(note_line)
                lines.append("")
            
            # Interpretation
            lines.append(f"   These markers suggest genetic connection to {name}.")
        else:
            lines.append("   No characteristic markers detected.")
            lines.append("")
            total = pop_data.get("total_markers", 0)
            if total > 0:
                lines.append(f"   (Checked {total} markers, none matched)")
        
        # Modern contribution
        contrib = metadata.get("modern_contribution", "")
        if contrib:
            lines.append(f"\n   Modern contribution: {contrib}")
        
        lines.append("")
    
    # Migration narrative
    lines.append("=" * 70)
    lines.append("📜 THE STORY OF EUROPEAN ANCESTRY")
    lines.append("=" * 70)
    lines.append("")
    lines.append(_get_migration_narrative(signals))
    lines.append("")
    
    # Methodology
    lines.append("=" * 70)
    lines.append("📖 METHODOLOGY & LIMITATIONS")
    lines.append("=" * 70)
    lines.append("")
    lines.append("This analysis is based on published ancient DNA studies:")
    lines.append("• Lazaridis et al. 2014 (PMID: 25230663)")
    lines.append("• Haak et al. 2015 (PMID: 25731166)")
    lines.append("• Mathieson et al. 2015 (PMID: 26595274)")
    lines.append("• Prüfer et al. 2014 (PMID: 24352235)")
    lines.append("")
    lines.append("LIMITATIONS:")
    lines.append("• Consumer arrays only cover a small fraction of ancestry-informative markers")
    lines.append("• Individual marker presence/absence doesn't equal ancestry percentage")
    lines.append("• All modern Europeans have ancestry from multiple ancient populations")
    lines.append("• This is educational content, not a definitive ancestry determination")
    lines.append("")
    
    return "\n".join(lines)


def _get_migration_narrative(signals: Dict[str, Any]) -> str:
    """Generate a narrative about European ancestry based on detected signals."""
    
    # Check what was detected
    whg_detected = len(signals.get("WHG", {}).get("markers_detected", [])) > 0
    anf_detected = len(signals.get("ANF", {}).get("markers_detected", [])) > 0
    yamnaya_detected = len(signals.get("Yamnaya", {}).get("markers_detected", [])) > 0
    neanderthal_detected = len(signals.get("Neanderthal", {}).get("markers_detected", [])) > 0
    
    narrative_parts = []
    
    # Introduction
    narrative_parts.append(
        "Modern European ancestry is a mix of THREE main ancient populations, "
        "arriving at different times over the past 40,000 years:"
    )
    
    # WHG
    narrative_parts.append("\n")
    narrative_parts.append(
        "1. MESOLITHIC HUNTER-GATHERERS (~10,000 BCE): The first modern humans to "
        "settle Europe after the Ice Age. Surprisingly, they had DARK SKIN and BLUE EYES - "
        "a combination not seen in modern populations. They lived by hunting and "
        "gathering for thousands of years."
    )
    if whg_detected:
        narrative_parts.append(
            "   → You carry markers associated with this ancient population."
        )
    
    # ANF
    narrative_parts.append("\n")
    narrative_parts.append(
        "2. NEOLITHIC FARMERS (~7,000 BCE): Migrants from Anatolia (modern Turkey) who "
        "brought agriculture to Europe. They had LIGHT SKIN and BROWN EYES. Within "
        "2,000 years, they replaced or absorbed most hunter-gatherers."
    )
    if anf_detected:
        narrative_parts.append(
            "   → You carry markers associated with these early farmers."
        )
    
    # Yamnaya
    narrative_parts.append("\n")
    narrative_parts.append(
        "3. BRONZE AGE STEPPE HERDERS (~3,000 BCE): The Yamnaya people from the "
        "Pontic steppe migrated massively into Europe, mixing with existing populations. "
        "They likely brought Indo-European languages and the horse. They contributed "
        "30-50% of Northern European ancestry."
    )
    if yamnaya_detected:
        narrative_parts.append(
            "   → You carry markers associated with Steppe ancestry."
        )
    
    # Neanderthal
    if neanderthal_detected:
        narrative_parts.append("\n")
        narrative_parts.append(
            "4. ARCHAIC ADMIXTURE: You also carry markers from Neanderthal "
            "introgression - DNA inherited from interbreeding between modern humans "
            "and Neanderthals ~50,000 years ago. This ~1-4% Neanderthal DNA affects "
            "immunity, skin, and other traits."
        )
    
    return "\n".join(narrative_parts)


def get_ancient_dna_json(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Get ancient DNA analysis data as structured JSON for dashboard.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Structured data for dashboard visualization
    """
    signals = detect_ancient_signals(genotypes)
    
    # Format for dashboard
    populations = []
    for pop_key in ["WHG", "EHG", "ANF", "Yamnaya", "Neanderthal", "Denisovan"]:
        pop_data = signals.get(pop_key, {})
        metadata = ANCIENT_POPULATIONS.get(pop_key, {})
        
        populations.append({
            "id": pop_key,
            "name": metadata.get("name", pop_key),
            "emoji": metadata.get("emoji", ""),
            "period": metadata.get("period", ""),
            "region": metadata.get("region", ""),
            "description": metadata.get("description", ""),
            "physical": metadata.get("physical", ""),
            "key_sites": metadata.get("key_sites", []),
            "modern_contribution": metadata.get("modern_contribution", ""),
            "pmid": metadata.get("pmid", []),
            "markers_detected": pop_data.get("markers_detected", []),
            "total_markers_checked": pop_data.get("total_markers", 0),
            "detection_count": pop_data.get("detection_count", 0),
        })
    
    return {
        "populations": populations,
        "narrative": _get_migration_narrative(signals),
        "timeline": [
            {"year": -40000, "event": "Modern humans arrive in Europe"},
            {"year": -15000, "event": "Ice Age ends, WHG populations expand"},
            {"year": -7000, "event": "Neolithic farmers arrive from Anatolia"},
            {"year": -3000, "event": "Yamnaya migrations from the Steppe"},
            {"year": -2000, "event": "Bronze Age - modern European genetics established"},
        ],
        "methodology": {
            "sources": [
                {"pmid": "25230663", "title": "Lazaridis et al. 2014 - Ancient human genomes"},
                {"pmid": "25731166", "title": "Haak et al. 2015 - Massive migration from the steppe"},
                {"pmid": "26595274", "title": "Mathieson et al. 2015 - Genome-wide patterns of selection"},
                {"pmid": "24352235", "title": "Prüfer et al. 2014 - Complete Neanderthal genome"},
            ],
            "limitations": [
                "Consumer arrays detect only a fraction of ancestry-informative markers",
                "Individual marker detection doesn't equal ancestry percentage",
                "All modern Europeans have mixed ancient ancestry",
                "This is educational, not a definitive ancestry determination"
            ]
        }
    }


# =============================================================================
# NEANDERTHAL-SPECIFIC ANALYSIS
# =============================================================================

def get_neanderthal_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Generate detailed Neanderthal introgression analysis.
    
    Args:
        genotypes: User's genotypes
        
    Returns:
        Neanderthal-specific analysis
    """
    ancient_data = get_ancient_data()
    
    neanderthal_markers = []
    total_checked = 0
    
    for rsid, marker_data in ancient_data.items():
        if rsid.startswith("_"):
            continue
        
        introgression = marker_data.get("introgression", {})
        if introgression.get("source") != "Neanderthal":
            continue
        
        total_checked += 1
        user_geno = genotypes.get(rsid)
        if not user_geno:
            continue
        
        derived = marker_data.get("derived_allele", "")
        has_neanderthal = derived and derived.upper() in user_geno.upper()
        
        if has_neanderthal:
            neanderthal_markers.append({
                "rsid": rsid,
                "gene": marker_data.get("gene", ""),
                "trait": marker_data.get("trait", ""),
                "your_genotype": user_geno,
                "neanderthal_allele": derived,
                "effect": introgression.get("effect", ""),
                "modern_frequency_EUR": introgression.get("modern_frequency_EUR", 0),
                "history": marker_data.get("history", "")
            })
    
    # Estimate based on detected markers (very rough)
    detection_rate = len(neanderthal_markers) / total_checked if total_checked > 0 else 0
    
    return {
        "markers_with_neanderthal_variants": neanderthal_markers,
        "total_neanderthal_markers_checked": total_checked,
        "detected_count": len(neanderthal_markers),
        "detection_rate": round(detection_rate * 100, 1),
        "context": {
            "typical_european": "~1.5-2.5% Neanderthal DNA",
            "typical_east_asian": "~1.5-2.5% Neanderthal DNA",
            "typical_african": "~0% Neanderthal DNA (except through recent admixture)",
            "note": "Consumer arrays can only detect a small fraction of Neanderthal variants"
        },
        "interpretation": _interpret_neanderthal(len(neanderthal_markers), total_checked),
        "pmid": ["24352235", "23128226", "29532074"]
    }


def _interpret_neanderthal(detected: int, total: int) -> str:
    """Generate interpretation of Neanderthal detection."""
    if total == 0:
        return "No Neanderthal markers available for analysis."
    
    if detected == 0:
        return (
            "No Neanderthal variants detected among the markers checked. "
            "This doesn't mean you have no Neanderthal ancestry - consumer arrays "
            "only test a tiny fraction of the ~40,000 Neanderthal variants in modern humans."
        )
    elif detected <= 2:
        return (
            f"Detected {detected} Neanderthal-derived variant(s). This is consistent "
            "with typical non-African ancestry (~1-4% Neanderthal DNA). These variants "
            "were inherited from interbreeding events ~50,000 years ago."
        )
    else:
        return (
            f"Detected {detected} Neanderthal-derived variants. You carry several "
            "introgressed regions from Neanderthals, which is normal for non-African populations. "
            "Some of these variants may have been positively selected (like immune genes)."
        )
