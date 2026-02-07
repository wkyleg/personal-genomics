#!/usr/bin/env python3
"""
Test statistical rigor against Kyle's DNA data.

This script demonstrates the new statistical features:
- Confidence intervals
- P-values
- Confidence levels
"""

import json
import sys
from pathlib import Path

# Add the skill directory to path
sys.path.insert(0, str(Path(__file__).parent))

from personal_genomics.statistics import (
    ConfidenceLevel, wilson_score_interval, 
    prs_percentile_ci, ancestry_similarity_stats
)
from markers.population_comparison import find_most_similar_populations
from markers.polygenic_scores import calculate_prs, PRS_CONDITIONS
from markers.ancient_ancestry import detect_ancient_signals
from markers.pharmacogenomics_stats import analyze_all_pharmacogenes
from markers.trait_stats import predict_all_traits


def load_kyle_genotypes(filepath: str = "/Users/bran/dna-analysis/raw_data.txt") -> dict:
    """Load Kyle's genotypes from raw file."""
    genotypes = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                rsid = parts[0]
                allele1 = parts[3]
                allele2 = parts[4]
                
                if rsid.startswith('rs'):
                    genotypes[rsid] = allele1 + allele2
    
    return genotypes


def main():
    print("=" * 70)
    print("STATISTICAL RIGOR TEST - KYLE'S DNA")
    print("=" * 70)
    print()
    
    # Load genotypes
    print("Loading genotypes...")
    genotypes = load_kyle_genotypes()
    print(f"Loaded {len(genotypes):,} SNPs")
    print()
    
    # Test 1: Population Comparison with Statistics
    print("-" * 70)
    print("1. POPULATION COMPARISON (with CIs and p-values)")
    print("-" * 70)
    
    similar_pops = find_most_similar_populations(genotypes, top_n=5)
    
    for pop in similar_pops:
        print(f"\n{pop['flag']} {pop['population_name']} ({pop['superpopulation']})")
        print(f"   Similarity: {pop['similarity']:.1f}% (95% CI: {pop['ci_lower']:.1f}-{pop['ci_upper']:.1f}%)")
        if pop.get('p_value'):
            print(f"   P-value vs random: {pop['p_value']:.4f}")
        print(f"   Based on {pop['n_markers']} markers")
        print(f"   Confidence: {pop['confidence']}")
    
    # Test 2: PRS with Confidence Intervals
    print("\n" + "-" * 70)
    print("2. POLYGENIC RISK SCORES (with CIs)")
    print("-" * 70)
    
    for condition in ["cad", "t2d", "alzheimer", "obesity"]:
        prs = calculate_prs(genotypes, condition)
        
        # Get display name from PRS_CONDITIONS
        condition_info = {k: v for k, v in PRS_CONDITIONS.items() 
                        if condition in k.lower() or k.lower().replace("_", "") == condition}
        display_name = list(condition_info.keys())[0].replace("_", " ").title() if condition_info else condition.upper()
        
        print(f"\n{display_name}:")
        print(f"   Percentile: {prs['percentile']} (95% CI: {prs['ci_lower']}-{prs['ci_upper']})")
        print(f"   Based on {prs['markers_found']} of {prs['markers_total']} markers")
        print(f"   Confidence: {prs['confidence']}")
        if prs.get('warnings'):
            for w in prs['warnings']:
                print(f"   ⚠️ {w}")
    
    # Test 3: Ancient Ancestry with Statistics
    print("\n" + "-" * 70)
    print("3. ANCIENT ANCESTRY SIGNALS (with CIs)")
    print("-" * 70)
    
    ancient = detect_ancient_signals(genotypes)
    
    for pop in ["WHG", "ANF", "Yamnaya", "Neanderthal"]:
        if pop in ancient:
            data = ancient[pop]
            print(f"\n{data.get('metadata', {}).get('emoji', '')} {data.get('metadata', {}).get('name', pop)}:")
            print(f"   Signal: {data.get('signal', 'N/A')}")
            print(f"   Derived ratio: {data.get('derived_ratio', 0):.3f} (95% CI: {data.get('ci_lower', 0):.3f}-{data.get('ci_upper', 0):.3f})")
            print(f"   Markers detected: {data.get('detection_count', 0)} of {data.get('n_markers', 0)}")
            print(f"   Confidence: {data.get('confidence', 'N/A')}")
    
    # Test 4: Pharmacogenomics with Confidence
    print("\n" + "-" * 70)
    print("4. PHARMACOGENOMICS (with diplotype confidence)")
    print("-" * 70)
    
    pgx = analyze_all_pharmacogenes(genotypes)
    
    for gene, result in pgx.items():
        r = result.to_dict()
        print(f"\n{gene}:")
        print(f"   Phenotype: {r['phenotype']}")
        print(f"   Diplotype: {r['diplotype']}")
        print(f"   Activity score: {r['activity_score']:.1f} ± {r['activity_uncertainty']:.2f}")
        print(f"   Confidence: {r['confidence']} ({r['diplotype_confidence']:.0%})")
        if r['warnings']:
            for w in r['warnings']:
                print(f"   ⚠️ {w}")
    
    # Test 5: Traits with Probability CIs
    print("\n" + "-" * 70)
    print("5. TRAIT PREDICTIONS (with probability CIs)")
    print("-" * 70)
    
    traits = predict_all_traits(genotypes)
    
    for rsid, pred in list(traits.items())[:6]:
        p = pred.to_dict()
        print(f"\n{p['trait']} ({p['gene']}):")
        print(f"   Probability: {p['probability']*100:.1f}% (95% CI: {p['ci_lower']*100:.1f}-{p['ci_upper']*100:.1f}%)")
        print(f"   Odds ratio: {p['odds_ratio']:.1f} (95% CI: {p['or_ci_lower']:.1f}-{p['or_ci_upper']:.1f})")
        print(f"   Genotype: {p['genotype']}")
        print(f"   Confidence: {p['confidence']}")
    
    print("\n" + "=" * 70)
    print("STATISTICAL RIGOR TEST COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
