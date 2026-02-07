#!/usr/bin/env python3
"""
Test script for v5.0 genomics expansion.
Tests against Kyle's DNA data.
"""

import sys
import json
from pathlib import Path

# Add markers module to path
sys.path.insert(0, str(Path(__file__).parent))

from markers.v5_integration import (
    V5_ALL_MARKERS,
    get_v5_marker_counts,
    parse_dna_file,
    generate_comprehensive_v5_report
)

def main():
    dna_path = Path.home() / "dna-analysis" / "raw_data.txt"
    
    if not dna_path.exists():
        print(f"❌ DNA file not found: {dna_path}")
        return 1
    
    print("🧬 Testing v5.0 Genomics Expansion")
    print("=" * 60)
    
    # Get marker counts
    counts = get_v5_marker_counts()
    print("\n📊 V5.0 Marker Counts by Category:")
    for category, count in counts.items():
        print(f"  {category}: {count}")
    
    # Parse DNA
    print(f"\n📁 Parsing DNA file: {dna_path}")
    genotypes = parse_dna_file(str(dna_path))
    print(f"  Total SNPs in file: {len(genotypes)}")
    
    # Check coverage
    v5_found = sum(1 for rs in V5_ALL_MARKERS if rs in genotypes)
    print(f"  V5 markers found: {v5_found}/{counts['total_unique']} ({100*v5_found/counts['total_unique']:.1f}%)")
    
    # Generate report
    print("\n🔬 Generating comprehensive report...")
    report = generate_comprehensive_v5_report(genotypes)
    
    # Print critical findings
    print("\n🚨 CRITICAL FINDINGS:")
    if report.get("critical_findings_summary"):
        for finding in report["critical_findings_summary"]:
            print(f"  {finding}")
    else:
        print("  ✅ No critical findings")
    
    # Print key results
    print("\n📋 KEY RESULTS:")
    
    # Pharmacogenomics
    pharma = report.get("pharmacogenomics", {})
    print(f"\n  Pharmacogenomics:")
    print(f"    Markers analyzed: {pharma.get('markers_analyzed', 0)}")
    if "metabolizer_statuses" in pharma:
        for gene, status in pharma["metabolizer_statuses"].items():
            print(f"    {gene}: {status.get('status', 'unknown')}")
    
    # Daily optimization
    daily = report.get("daily_optimization", {})
    chrono = daily.get("chronotype", {})
    print(f"\n  Daily Optimization:")
    print(f"    Chronotype: {chrono.get('type', 'unknown')}")
    
    caffeine = daily.get("caffeine", {})
    print(f"    Caffeine metabolism: {caffeine.get('metabolism', 'unknown')}")
    print(f"    Caffeine cutoff: {caffeine.get('caffeine_cutoff_time', 'unknown')}")
    
    # Athletic
    athletic = report.get("athletic", {})
    profile = athletic.get("athletic_profile", {})
    print(f"\n  Athletic Profile:")
    print(f"    Type: {profile.get('profile', 'unknown')}")
    print(f"    ACTN3: {profile.get('actn3_status', 'unknown')}")
    
    # Nutrition
    nutrition = report.get("nutrition", {})
    print(f"\n  Nutrition:")
    print(f"    APOE: {nutrition.get('apoe_genotype', 'unknown')}")
    diet = nutrition.get("diet_recommendations", {})
    print(f"    Diet type: {diet.get('type', 'unknown')}")
    
    # Longevity
    longevity = report.get("longevity", {})
    print(f"\n  Longevity:")
    print(f"    Score: {longevity.get('longevity_score', 'N/A')}")
    print(f"    Potential: {longevity.get('potential', 'unknown')}")
    
    # Mental health
    mental = report.get("mental_health", {})
    stress = mental.get("stress_phenotype", {})
    print(f"\n  Mental Health/Stress:")
    print(f"    COMT type: {stress.get('stress_type', 'unknown')}")
    
    # Quirky traits
    quirky = report.get("quirky_traits", {})
    print(f"\n  Quirky Traits:")
    for trait in quirky.get("quirky_traits", [])[:5]:
        print(f"    {trait.get('trait')}: {trait.get('result')}")
    
    # Medical special
    medical = report.get("medical_special", {})
    print(f"\n  Medical Special:")
    celiac = medical.get("celiac_hla", {})
    print(f"    Celiac HLA risk: {celiac.get('risk_level', 'unknown')}")
    blood = medical.get("blood_type", {})
    print(f"    Inferred blood type: {blood.get('inferred_abo', 'unknown')}")
    
    infection = medical.get("infection_resistance", {})
    if infection.get("resistances"):
        print(f"    Infection resistances:")
        for r in infection["resistances"]:
            print(f"      - {r['pathogen']}: {r['status']}")
    
    # Save full report
    report_path = Path(__file__).parent / "kyle_v5_report.json"
    with open(report_path, 'w') as f:
        # Convert enums to strings for JSON serialization
        def convert(obj):
            if hasattr(obj, 'value'):
                return obj.value
            raise TypeError(f"Cannot serialize {type(obj)}")
        json.dump(report, f, indent=2, default=convert)
    print(f"\n💾 Full report saved to: {report_path}")
    
    print("\n✅ Test completed successfully!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
