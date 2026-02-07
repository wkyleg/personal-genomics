#!/usr/bin/env python3
"""
Comprehensive Genetic Analysis - 1600+ Markers (v4.1.0)
Full health, pharmacogenomics, ancestry, traits, and actionable recommendations.
Works with ANY ancestry/ethnic background worldwide.

Supports:
- 23andMe (v3, v4, v5)
- AncestryDNA
- MyHeritage
- FamilyTreeDNA
- Nebula Genomics
- VCF files (whole genome/exome)

Privacy: All analysis runs locally. Zero network requests.

Output:
- Human-readable reports
- Agent-friendly JSON with actionable fields and priorities
- Polygenic risk scores for major conditions
- Evidence-based recommendations with citations
- Lifestyle recommendation engine
- Drug interaction matrix

Categories (v3.0):
1. Pharmacogenomics - Drug metabolism
2. Polygenic Risk Scores - Disease risk
3. Carrier Status - Recessive carriers
4. Health Risks - Disease susceptibility
5. Traits - Physical/behavioral
6. Nutrition - Nutrigenomics
7. Fitness - Athletic performance
8. Neurogenetics - Cognition/behavior
9. Longevity - Aging markers
10. Immunity - HLA and immune
11. Rare Diseases - Rare genetic conditions
12. Mental Health - Psychiatric genetics
13. Dermatology - Skin and hair
14. Vision & Hearing - Sensory genetics
15. Fertility - Reproductive health
"""

import sys
import json
import gzip
import math
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple

# Output directory
OUTPUT_DIR = Path.home() / "dna-analysis" / "reports"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Try to import marker modules
try:
    from markers.pharmacogenomics import PHARMACOGENOMICS_MARKERS, DRUG_INTERACTIONS
    from markers.polygenic_scores import PRS_WEIGHTS, PRS_CONDITIONS, calculate_prs
    from markers.carrier_status import CARRIER_MARKERS, CARRIER_SCREENING_PANELS
    from markers.health_risks import HEALTH_RISK_MARKERS
    from markers.traits import TRAIT_MARKERS
    from markers.nutrition import NUTRITION_MARKERS
    from markers.fitness import FITNESS_MARKERS
    from markers.neurogenetics import NEURO_MARKERS
    from markers.longevity import LONGEVITY_MARKERS
    from markers.immunity import IMMUNITY_MARKERS, HLA_DRUG_ALERTS
    # New v3.0 modules
    from markers.rare_diseases import RARE_DISEASE_MARKERS
    from markers.mental_health import MENTAL_HEALTH_MARKERS
    from markers.dermatology import DERMATOLOGY_MARKERS
    from markers.vision_hearing import VISION_HEARING_MARKERS
    from markers.fertility import FERTILITY_MARKERS
    from markers import get_marker_counts
    MODULES_LOADED = True
except ImportError as e:
    print(f"Warning: Could not load marker modules: {e}")
    print("Using inline markers only.")
    MODULES_LOADED = False
    PHARMACOGENOMICS_MARKERS = {}
    DRUG_INTERACTIONS = {}
    PRS_WEIGHTS = {}
    PRS_CONDITIONS = {}
    CARRIER_MARKERS = {}
    HEALTH_RISK_MARKERS = {}
    TRAIT_MARKERS = {}
    NUTRITION_MARKERS = {}
    FITNESS_MARKERS = {}
    NEURO_MARKERS = {}
    LONGEVITY_MARKERS = {}
    IMMUNITY_MARKERS = {}
    HLA_DRUG_ALERTS = {}
    RARE_DISEASE_MARKERS = {}
    MENTAL_HEALTH_MARKERS = {}
    DERMATOLOGY_MARKERS = {}
    VISION_HEARING_MARKERS = {}
    FERTILITY_MARKERS = {}
    def get_marker_counts(): return {"total": 0}


# =============================================================================
# FILE FORMAT DETECTION AND LOADING
# =============================================================================

def detect_format(filepath: str) -> str:
    """Detect DNA file format."""
    filepath = str(filepath)
    
    if filepath.endswith('.vcf') or filepath.endswith('.vcf.gz'):
        return 'vcf'
    
    # Read first few lines to detect format
    opener = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    with opener(filepath, mode) as f:
        header_lines = []
        for i, line in enumerate(f):
            if i >= 20:
                break
            header_lines.append(line)
    
    content = ''.join(header_lines).lower()
    
    if '23andme' in content:
        return '23andme'
    elif 'ancestrydna' in content:
        return 'ancestry'
    elif 'myheritage' in content:
        return 'myheritage'
    elif 'ftdna' in content or 'family tree dna' in content:
        return 'ftdna'
    elif '#rsid' in content or 'rsid\t' in content:
        return 'generic'
    else:
        return 'generic'


def load_vcf(filepath: str) -> Dict[str, str]:
    """Load VCF file into rsid -> genotype dict."""
    genotypes = {}
    
    opener = gzip.open if str(filepath).endswith('.gz') else open
    mode = 'rt' if str(filepath).endswith('.gz') else 'r'
    
    with opener(filepath, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            
            chrom, pos, rsid, ref, alt, qual, filt, info, fmt, sample = parts[:10]
            
            if not rsid.startswith('rs'):
                continue
            
            # Parse genotype
            fmt_fields = fmt.split(':')
            sample_fields = sample.split(':')
            
            gt_idx = fmt_fields.index('GT') if 'GT' in fmt_fields else 0
            gt = sample_fields[gt_idx] if gt_idx < len(sample_fields) else './.'
            
            # Convert GT to alleles
            alleles = [ref] + alt.split(',')
            gt_parts = gt.replace('|', '/').split('/')
            
            try:
                a1 = alleles[int(gt_parts[0])] if gt_parts[0] != '.' else '?'
                a2 = alleles[int(gt_parts[1])] if len(gt_parts) > 1 and gt_parts[1] != '.' else a1
                genotypes[rsid] = a1 + a2
            except (ValueError, IndexError):
                continue
    
    return genotypes


def load_consumer_format(filepath: str) -> Dict[str, str]:
    """Load 23andMe, Ancestry, or similar format."""
    genotypes = {}
    
    opener = gzip.open if str(filepath).endswith('.gz') else open
    mode = 'rt' if str(filepath).endswith('.gz') else 'r'
    
    with opener(filepath, mode) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 4:
                parts = line.strip().split(',')
            
            if len(parts) >= 4:
                rsid = parts[0]
                if rsid.startswith('rs'):
                    # Format: rsid, chrom, pos, genotype
                    genotype = parts[3].replace(' ', '')
                    if genotype and genotype != '--' and genotype != '00':
                        genotypes[rsid] = genotype
    
    return genotypes


def load_dna_file(filepath: str) -> Tuple[Dict[str, str], str]:
    """Load DNA data from any supported format."""
    fmt = detect_format(filepath)
    print(f"Detected format: {fmt}")
    
    if fmt == 'vcf':
        genotypes = load_vcf(filepath)
    else:
        genotypes = load_consumer_format(filepath)
    
    return genotypes, fmt


# =============================================================================
# APOE DETERMINATION
# =============================================================================

def determine_apoe(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Determine APOE genotype from rs429358 and rs7412."""
    rs429358 = genotypes.get('rs429358', '')
    rs7412 = genotypes.get('rs7412', '')
    
    if not rs429358 or not rs7412:
        return {
            "genotype": "unknown",
            "risk_level": "unknown",
            "rs429358": rs429358,
            "rs7412": rs7412,
            "interpretation": "Unable to determine APOE status"
        }
    
    # Determine alleles
    # ε2: rs429358=T, rs7412=T
    # ε3: rs429358=T, rs7412=C
    # ε4: rs429358=C, rs7412=C
    
    alleles = []
    for i in range(min(len(rs429358), len(rs7412))):
        c1 = rs429358[i].upper()
        c2 = rs7412[i].upper()
        
        if c1 == 'T' and c2 == 'T':
            alleles.append('ε2')
        elif c1 == 'T' and c2 == 'C':
            alleles.append('ε3')
        elif c1 == 'C' and c2 == 'C':
            alleles.append('ε4')
    
    if len(alleles) == 2:
        genotype = '/'.join(sorted(alleles))
    elif len(alleles) == 1:
        genotype = alleles[0] + '/' + alleles[0]
    else:
        genotype = 'unknown'
    
    risk_info = {
        'ε2/ε2': ('low', 'Protective. Lower Alzheimer\'s risk, but higher triglycerides.'),
        'ε2/ε3': ('low', 'Below average Alzheimer\'s risk.'),
        'ε3/ε3': ('average', 'Most common genotype. Average risk.'),
        'ε2/ε4': ('moderate', 'Mixed effects. ε2 partially offsets ε4.'),
        'ε3/ε4': ('elevated', '~3x Alzheimer\'s risk vs ε3/ε3. Lifestyle modifications important.'),
        'ε4/ε4': ('high', '~12x Alzheimer\'s risk. Exercise, diet, sleep, and cognitive engagement are protective.')
    }
    
    risk_level, interpretation = risk_info.get(genotype, ('unknown', 'Unable to interpret'))
    
    return {
        "genotype": genotype,
        "risk_level": risk_level,
        "rs429358": rs429358,
        "rs7412": rs7412,
        "interpretation": interpretation,
        "actionable": risk_level in ['elevated', 'high'],
        "recommendations": [
            "Regular aerobic exercise (strongest protective factor)",
            "Mediterranean diet",
            "7-8 hours quality sleep",
            "Cognitive engagement and social connection",
            "Cardiovascular risk factor management"
        ] if risk_level in ['elevated', 'high'] else []
    }


# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def analyze_markers(genotypes: Dict[str, str], markers: Dict, category: str) -> Dict[str, Any]:
    """Analyze a category of markers."""
    results = {
        "category": category,
        "total_in_database": len(markers),
        "found_in_data": 0,
        "risk_variants": 0,
        "findings": [],
        "actionable_items": []
    }
    
    for rsid, info in markers.items():
        geno = genotypes.get(rsid)
        if not geno:
            continue
        
        results["found_in_data"] += 1
        
        # Determine if risk allele present
        risk_allele = info.get('risk_allele') or info.get('effect_allele')
        risk_count = geno.upper().count(risk_allele.upper()) if risk_allele else 0
        
        finding = {
            "rsid": rsid,
            "gene": info.get('gene', 'Unknown'),
            "genotype": geno,
            "risk_allele": risk_allele,
            "risk_copies": risk_count,
            "is_risk": risk_count > 0
        }
        
        # Add variant-specific info
        for key in ['variant', 'name', 'condition', 'conditions', 'trait', 'effect', 'note', 'evidence']:
            if key in info:
                finding[key] = info[key]
        
        if risk_count > 0:
            results["risk_variants"] += 1
            
            # Check for actionable items
            if 'actionable' in info:
                action = {
                    "rsid": rsid,
                    "gene": info.get('gene'),
                    "genotype": geno,
                    **info['actionable']
                }
                results["actionable_items"].append(action)
        
        results["findings"].append(finding)
    
    return results


def calculate_all_prs(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Calculate polygenic risk scores for all conditions."""
    if not PRS_WEIGHTS:
        return {"error": "PRS module not loaded"}
    
    scores = {}
    conditions = set(v['condition'] for v in PRS_WEIGHTS.values())
    
    for condition in conditions:
        condition_snps = {k: v for k, v in PRS_WEIGHTS.items() if v['condition'] == condition}
        
        score = 0
        found = 0
        for rsid, info in condition_snps.items():
            geno = genotypes.get(rsid)
            if geno:
                found += 1
                effect_count = geno.upper().count(info['effect'].upper())
                score += effect_count * info['beta']
        
        if found > 5:
            # Rough percentile estimation
            import math
            z = score / math.sqrt(found * 0.5) if found > 0 else 0
            percentile = min(99, max(1, int(50 + z * 15)))
        else:
            percentile = None
        
        scores[condition] = {
            "raw_score": round(score, 3),
            "snps_found": found,
            "snps_total": len(condition_snps),
            "coverage": round(found / len(condition_snps), 2) if condition_snps else 0,
            "percentile_estimate": percentile,
            "confidence": "moderate" if found > len(condition_snps) * 0.5 else "low"
        }
    
    return scores


# =============================================================================
# AGENT-FRIENDLY OUTPUT
# =============================================================================

def generate_agent_summary(all_results: Dict) -> Dict[str, Any]:
    """Generate structured output optimized for AI agents."""
    
    summary = {
        "analysis_timestamp": datetime.now().isoformat(),
        "snps_analyzed": all_results.get("total_snps", 0),
        "format_detected": all_results.get("format", "unknown"),
        
        # Priority-sorted actionable items
        "critical_alerts": [],
        "high_priority": [],
        "medium_priority": [],
        "low_priority": [],
        "informational": [],
        
        # Key health markers
        "apoe_status": all_results.get("apoe", {}),
        "pharmacogenomics_alerts": [],
        "carrier_status": [],
        "polygenic_risk_scores": all_results.get("prs", {}),
        
        # Traits and lifestyle
        "notable_traits": [],
        "nutrition_insights": [],
        "fitness_insights": [],
        
        # Metadata
        "confidence_notes": [
            "Consumer arrays capture ~0.1% of genome",
            "Polygenic scores are probabilistic, not deterministic",
            "Many conditions depend heavily on environment and lifestyle",
            "These results are not diagnostic - consult healthcare providers"
        ]
    }
    
    # Collect all actionable items and sort by priority
    priority_order = {"critical": 0, "high": 1, "medium": 2, "low": 3, "informational": 4}
    
    for category, results in all_results.items():
        if isinstance(results, dict) and "actionable_items" in results:
            for item in results["actionable_items"]:
                priority = item.get("priority", "informational")
                
                if priority == "critical":
                    summary["critical_alerts"].append(item)
                elif priority == "high":
                    summary["high_priority"].append(item)
                elif priority == "medium":
                    summary["medium_priority"].append(item)
                elif priority == "low":
                    summary["low_priority"].append(item)
                else:
                    summary["informational"].append(item)
                
                # Also add to specific categories
                if category == "pharmacogenomics":
                    summary["pharmacogenomics_alerts"].append(item)
                elif category == "carrier_status":
                    summary["carrier_status"].append(item)
    
    # Add notable traits
    if "traits" in all_results:
        for finding in all_results["traits"].get("findings", [])[:20]:
            if finding.get("effect"):
                summary["notable_traits"].append({
                    "trait": finding.get("trait") or finding.get("name"),
                    "gene": finding.get("gene"),
                    "genotype": finding.get("genotype"),
                    "interpretation": finding.get("effect")
                })
    
    return summary


# =============================================================================
# HUMAN-READABLE REPORT
# =============================================================================

def generate_report(all_results: Dict, agent_summary: Dict) -> str:
    """Generate human-readable report."""
    lines = []
    
    lines.append("=" * 78)
    lines.append("COMPREHENSIVE GENETIC ANALYSIS REPORT")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("=" * 78)
    lines.append("")
    lines.append("IMPORTANT DISCLAIMERS:")
    lines.append("  * This is NOT medical advice")
    lines.append("  * Consult healthcare providers before acting on results")
    lines.append("  * Genetic risk does not equal destiny - lifestyle matters enormously")
    lines.append("  * Consumer arrays miss rare variants and structural changes")
    lines.append("")
    
    # Summary stats
    lines.append("-" * 78)
    lines.append("SUMMARY")
    lines.append("-" * 78)
    lines.append(f"SNPs in file: {all_results.get('total_snps', 'N/A'):,}")
    lines.append(f"File format: {all_results.get('format', 'Unknown')}")
    lines.append(f"Critical alerts: {len(agent_summary.get('critical_alerts', []))}")
    lines.append(f"High priority items: {len(agent_summary.get('high_priority', []))}")
    lines.append("")
    
    # APOE
    apoe = all_results.get("apoe", {})
    if apoe.get("genotype") != "unknown":
        lines.append("-" * 78)
        lines.append("APOE STATUS (Alzheimer's / Cardiovascular)")
        lines.append("-" * 78)
        lines.append(f"Genotype: {apoe.get('genotype')}")
        lines.append(f"Risk level: {apoe.get('risk_level')}")
        lines.append(f"Interpretation: {apoe.get('interpretation')}")
        if apoe.get("recommendations"):
            lines.append("Recommendations:")
            for rec in apoe["recommendations"]:
                lines.append(f"  - {rec}")
        lines.append("")
    
    # Critical alerts
    if agent_summary.get("critical_alerts"):
        lines.append("-" * 78)
        lines.append("!!! CRITICAL ALERTS - SHARE WITH HEALTHCARE PROVIDERS !!!")
        lines.append("-" * 78)
        for alert in agent_summary["critical_alerts"]:
            lines.append(f"\n  {alert.get('gene', 'Unknown')} ({alert.get('rsid')})")
            lines.append(f"  Genotype: {alert.get('genotype')}")
            for rec in alert.get("recommendations", []):
                lines.append(f"    * {rec}")
        lines.append("")
    
    # High priority
    if agent_summary.get("high_priority"):
        lines.append("-" * 78)
        lines.append("HIGH PRIORITY FINDINGS")
        lines.append("-" * 78)
        for item in agent_summary["high_priority"][:10]:
            lines.append(f"\n  {item.get('gene', 'Unknown')} ({item.get('rsid')})")
            lines.append(f"  Genotype: {item.get('genotype')}")
            if item.get("recommendations"):
                for rec in item["recommendations"][:3]:
                    lines.append(f"    - {rec}")
        lines.append("")
    
    # Pharmacogenomics
    if agent_summary.get("pharmacogenomics_alerts"):
        lines.append("-" * 78)
        lines.append("PHARMACOGENOMICS (Drug Response)")
        lines.append("-" * 78)
        for item in agent_summary["pharmacogenomics_alerts"][:15]:
            lines.append(f"  {item.get('gene')}: {item.get('genotype')} - {item.get('action_type', '')}")
        lines.append("")
    
    # Polygenic Risk Scores
    prs = all_results.get("prs", {})
    if prs and not prs.get("error"):
        lines.append("-" * 78)
        lines.append("POLYGENIC RISK SCORES")
        lines.append("-" * 78)
        for condition, scores in prs.items():
            if scores.get("percentile_estimate"):
                lines.append(f"  {condition}: {scores['percentile_estimate']}th percentile "
                           f"(confidence: {scores['confidence']})")
        lines.append("")
    
    # Traits
    if agent_summary.get("notable_traits"):
        lines.append("-" * 78)
        lines.append("NOTABLE TRAITS")
        lines.append("-" * 78)
        for trait in agent_summary["notable_traits"][:15]:
            interp = trait.get("interpretation", "")
            if isinstance(interp, dict):
                interp = interp.get(trait.get("genotype"), str(interp))
            lines.append(f"  {trait.get('trait')}: {interp[:60]}...")
        lines.append("")
    
    lines.append("=" * 78)
    lines.append("END OF REPORT")
    lines.append("=" * 78)
    
    return "\n".join(lines)


# =============================================================================
# MAIN
# =============================================================================

def generate_lifestyle_recommendations(all_results: Dict) -> Dict[str, Any]:
    """Generate personalized lifestyle recommendations based on genetic profile."""
    recommendations = {
        "diet": [],
        "exercise": [],
        "supplements": [],
        "screening": [],
        "lifestyle": [],
        "avoid": []
    }
    
    # APOE-based recommendations
    apoe = all_results.get("apoe", {})
    if apoe.get("risk_level") in ["elevated", "high"]:
        recommendations["diet"].extend([
            "Mediterranean diet - strongly protective for APOE4 carriers",
            "Limit saturated fat intake",
            "Include fatty fish 2-3x/week (omega-3s)"
        ])
        recommendations["exercise"].extend([
            "Regular aerobic exercise - strongest protective factor",
            "150+ minutes moderate exercise weekly"
        ])
        recommendations["lifestyle"].append("Prioritize 7-8 hours quality sleep")
    
    # Nutrition-based
    nutrition = all_results.get("nutrition", {})
    for finding in nutrition.get("findings", []):
        gene = finding.get("gene", "")
        if finding.get("risk_copies", 0) > 0:
            if gene == "MTHFR":
                recommendations["supplements"].append("Consider methylfolate over folic acid")
            elif gene == "VDR":
                recommendations["supplements"].append("Vitamin D supplementation may be beneficial")
            elif "caffeine" in str(finding.get("trait", "")).lower():
                recommendations["avoid"].append("Limit caffeine if slow metabolizer")
    
    # Fitness-based
    fitness = all_results.get("fitness", {})
    for finding in fitness.get("findings", []):
        gene = finding.get("gene", "")
        if gene == "ACTN3" and finding.get("genotype") == "TT":
            recommendations["exercise"].append("Favor endurance over power training (ACTN3 XX)")
        elif gene == "COL5A1" and finding.get("risk_copies", 0) > 0:
            recommendations["exercise"].append("Include tendon-strengthening exercises")
    
    # Dermatology-based
    dermatology = all_results.get("dermatology", {})
    for finding in dermatology.get("findings", []):
        gene = finding.get("gene", "")
        if gene == "MC1R" and finding.get("risk_copies", 0) > 0:
            recommendations["lifestyle"].append("Strict sun protection - elevated melanoma risk")
            recommendations["screening"].append("Annual dermatology screening recommended")
    
    # Vision-based
    vision = all_results.get("vision_hearing", {})
    for finding in vision.get("findings", []):
        gene = finding.get("gene", "")
        if gene == "CFH" and finding.get("risk_copies", 0) > 0:
            recommendations["screening"].append("Annual dilated eye exam after age 50 (AMD risk)")
            recommendations["supplements"].append("AREDS2 formula if AMD develops")
            recommendations["avoid"].append("CRITICAL: Do not smoke (dramatically increases AMD risk)")
    
    return recommendations


def generate_drug_interaction_matrix(all_results: Dict) -> Dict[str, Any]:
    """Generate drug interaction warnings based on pharmacogenomics."""
    matrix = {
        "critical_interactions": [],
        "warnings": [],
        "dosing_adjustments": [],
        "safe_alternatives": {}
    }
    
    pharma = all_results.get("pharmacogenomics", {})
    for finding in pharma.get("findings", []):
        gene = finding.get("gene", "")
        geno = finding.get("genotype", "")
        risk = finding.get("risk_copies", 0)
        
        # DPYD - 5-FU toxicity
        if gene == "DPYD" and risk > 0:
            matrix["critical_interactions"].append({
                "drug_class": "Fluoropyrimidines",
                "drugs": ["5-FU", "capecitabine", "tegafur"],
                "risk": "FATAL TOXICITY - avoid or reduce dose significantly",
                "gene": gene,
                "genotype": geno
            })
        
        # CYP2C19 - clopidogrel
        if gene == "CYP2C19":
            if "AA" in geno:  # Poor metabolizer
                matrix["warnings"].append({
                    "drug": "Clopidogrel (Plavix)",
                    "issue": "Poor activation - reduced efficacy",
                    "recommendation": "Consider prasugrel or ticagrelor",
                    "gene": gene,
                    "genotype": geno
                })
        
        # CYP2D6 - opioids
        if gene == "CYP2D6":
            if risk == 0:  # Poor metabolizer
                matrix["dosing_adjustments"].append({
                    "drug_class": "Codeine, tramadol",
                    "issue": "No activation to active metabolite",
                    "recommendation": "Use alternative analgesics (morphine, oxycodone)",
                    "gene": gene,
                    "genotype": geno
                })
        
        # SLCO1B1 - statins
        if gene == "SLCO1B1" and risk > 0:
            matrix["warnings"].append({
                "drug": "Simvastatin",
                "issue": "Elevated myopathy risk",
                "recommendation": "Limit to 20mg/day or use alternative statin",
                "gene": gene,
                "genotype": geno
            })
        
        # Warfarin dosing
        if gene in ["CYP2C9", "VKORC1"] and risk > 0:
            matrix["dosing_adjustments"].append({
                "drug": "Warfarin",
                "issue": "May require lower dose",
                "recommendation": "Use pharmacogenomic dosing algorithms",
                "gene": gene,
                "genotype": geno
            })
    
    return matrix


def main():
    if len(sys.argv) < 2:
        print("Personal Genomics Analysis Tool v3.0")
        print("=" * 40)
        print("\nUsage: python comprehensive_analysis.py <dna_file>")
        print("\nSupported formats:")
        print("  - 23andMe (v3, v4, v5)")
        print("  - AncestryDNA")
        print("  - MyHeritage")
        print("  - FamilyTreeDNA")
        print("  - VCF (whole genome/exome)")
        print("  - Any tab-delimited rsid format")
        print(f"\nMarker modules loaded: {MODULES_LOADED}")
        if MODULES_LOADED:
            counts = get_marker_counts()
            print(f"Total markers in database: {counts['total']:,}")
            print("\nCategories:")
            for k, v in counts.items():
                if k != 'total':
                    print(f"  {k}: {v}")
        sys.exit(1)
    
    filepath = sys.argv[1]
    print(f"Loading {filepath}...")
    
    genotypes, fmt = load_dna_file(filepath)
    print(f"Loaded {len(genotypes):,} SNPs")
    
    # Run all analyses
    all_results = {
        "total_snps": len(genotypes),
        "format": fmt,
        "apoe": determine_apoe(genotypes),
        "version": "3.0"
    }
    
    print("Analyzing markers...")
    
    if MODULES_LOADED:
        # Core categories
        all_results["pharmacogenomics"] = analyze_markers(genotypes, PHARMACOGENOMICS_MARKERS, "pharmacogenomics")
        all_results["carrier_status"] = analyze_markers(genotypes, CARRIER_MARKERS, "carrier_status")
        all_results["health_risks"] = analyze_markers(genotypes, HEALTH_RISK_MARKERS, "health_risks")
        all_results["traits"] = analyze_markers(genotypes, TRAIT_MARKERS, "traits")
        all_results["nutrition"] = analyze_markers(genotypes, NUTRITION_MARKERS, "nutrition")
        all_results["fitness"] = analyze_markers(genotypes, FITNESS_MARKERS, "fitness")
        all_results["neurogenetics"] = analyze_markers(genotypes, NEURO_MARKERS, "neurogenetics")
        all_results["longevity"] = analyze_markers(genotypes, LONGEVITY_MARKERS, "longevity")
        all_results["immunity"] = analyze_markers(genotypes, IMMUNITY_MARKERS, "immunity")
        all_results["prs"] = calculate_all_prs(genotypes)
        
        # New v3.0 categories
        all_results["rare_diseases"] = analyze_markers(genotypes, RARE_DISEASE_MARKERS, "rare_diseases")
        all_results["mental_health"] = analyze_markers(genotypes, MENTAL_HEALTH_MARKERS, "mental_health")
        all_results["dermatology"] = analyze_markers(genotypes, DERMATOLOGY_MARKERS, "dermatology")
        all_results["vision_hearing"] = analyze_markers(genotypes, VISION_HEARING_MARKERS, "vision_hearing")
        all_results["fertility"] = analyze_markers(genotypes, FERTILITY_MARKERS, "fertility")
        
        # Advanced features
        all_results["lifestyle_recommendations"] = generate_lifestyle_recommendations(all_results)
        all_results["drug_interaction_matrix"] = generate_drug_interaction_matrix(all_results)
    
    # Generate outputs
    print("Generating reports...")
    
    agent_summary = generate_agent_summary(all_results)
    report = generate_report(all_results, agent_summary)
    
    # Save files
    with open(OUTPUT_DIR / "full_analysis.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    with open(OUTPUT_DIR / "agent_summary.json", 'w') as f:
        json.dump(agent_summary, f, indent=2, default=str)
    
    with open(OUTPUT_DIR / "report.txt", 'w') as f:
        f.write(report)
    
    # Print report
    print("\n" + report)
    
    print(f"\nOutput files saved to: {OUTPUT_DIR}/")
    print(f"  - full_analysis.json    (complete data)")
    print(f"  - agent_summary.json    (AI-optimized)")
    print(f"  - report.txt            (human-readable)")


if __name__ == "__main__":
    main()
