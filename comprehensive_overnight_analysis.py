#!/usr/bin/env python3
"""
Comprehensive Overnight DNA Analysis

This script performs a complete analysis using all 9 reference datasets:
1. 1000 Genomes - Population comparison
2. HGDP - Fine-scale ancestry
3. SGDP - Global diversity
4. Ancient DNA - Ancestral signals
5. gnomAD - Allele frequencies
6. ClinVar - Clinical variants
7. PharmGKB - Pharmacogenomics
8. PGS Catalog - Polygenic risk scores
9. GWAS Catalog - Trait associations

Author: OpenClaw AI
Date: 2026-02-07
"""

from __future__ import annotations

import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TypedDict, Union
from dataclasses import dataclass, field, asdict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(
            Path.home() / "dna-analysis" / "reports" / "overnight_2026-02-07" / "analysis.log"
        )
    ]
)
logger = logging.getLogger(__name__)

# Add the skill to path
SKILL_PATH = Path.home() / ".openclaw" / "workspace" / "skills" / "personal-genomics"
sys.path.insert(0, str(SKILL_PATH))
sys.path.insert(0, str(SKILL_PATH / "personal_genomics"))

from datasets import (
    ThousandGenomes,
    HGDP,
    SGDPDataset,
    AncientDNADataset,
    GnomAD,
    ClinVar,
    PharmGKB,
    PGSCatalog,
    GWASCatalog,
    download_all_datasets,
    get_dataset_status,
    THOUSAND_GENOMES_POPULATIONS,
    SUPERPOPULATIONS,
    HGDP_POPULATIONS,
    HGDP_REGIONS,
)


# =============================================================================
# TYPE DEFINITIONS
# =============================================================================

class PopulationResult(TypedDict):
    """Population similarity result."""
    population: str
    name: str
    region: str
    similarity: float


class AncientSignal(TypedDict):
    """Ancient ancestry signal."""
    population: str
    signal_strength: str
    markers_matched: int
    total_markers: int
    traits: List[str]


class ClinicalVariant(TypedDict):
    """Clinical variant finding."""
    rsid: str
    gene: str
    significance: str
    conditions: List[str]
    genotype: str
    pmids: List[str]


class PharmacogenomicResult(TypedDict):
    """Pharmacogenomic finding."""
    gene: str
    phenotype: str
    activity_score: float
    drugs_affected: List[str]
    recommendation: str


class PRSResult(TypedDict):
    """Polygenic risk score result."""
    trait: str
    raw_score: float
    percentile: float
    risk_category: str
    interpretation: str
    pmid: str


class GWASFinding(TypedDict):
    """GWAS trait association."""
    rsid: str
    trait: str
    risk_allele: str
    user_genotype: str
    effect_direction: str
    odds_ratio: Optional[float]
    pmid: str


@dataclass
class ComprehensiveAnalysisResult:
    """Complete analysis results."""
    # Metadata
    analysis_date: str = ""
    dna_file: str = ""
    total_snps: int = 0
    snps_analyzed: int = 0
    
    # Population analysis
    top_1kg_populations: List[PopulationResult] = field(default_factory=list)
    superpopulation_breakdown: Dict[str, float] = field(default_factory=dict)
    top_hgdp_populations: List[PopulationResult] = field(default_factory=list)
    hgdp_region_breakdown: Dict[str, float] = field(default_factory=dict)
    top_sgdp_populations: List[PopulationResult] = field(default_factory=list)
    sgdp_region_breakdown: Dict[str, float] = field(default_factory=dict)
    
    # Ancient ancestry
    ancient_signals: List[AncientSignal] = field(default_factory=list)
    neanderthal_markers: int = 0
    neanderthal_total: int = 0
    
    # Clinical
    clinvar_findings: List[ClinicalVariant] = field(default_factory=list)
    pathogenic_count: int = 0
    vus_count: int = 0
    
    # Pharmacogenomics
    pharmacogenomics: List[PharmacogenomicResult] = field(default_factory=list)
    actionable_pgx: int = 0
    
    # PRS
    prs_results: List[PRSResult] = field(default_factory=list)
    high_risk_conditions: List[str] = field(default_factory=list)
    
    # GWAS
    gwas_findings: List[GWASFinding] = field(default_factory=list)
    notable_traits: List[str] = field(default_factory=list)


# =============================================================================
# DNA FILE LOADING
# =============================================================================

def load_dna_file(filepath: str) -> Dict[str, str]:
    """
    Load genotype data from DNA file.
    
    Supports AncestryDNA, 23andMe, and similar tab-delimited formats.
    
    Args:
        filepath: Path to the DNA data file
        
    Returns:
        Dict mapping rsID to genotype (e.g., {"rs1426654": "AA"})
    """
    logger.info(f"Loading DNA file: {filepath}")
    genotypes: Dict[str, str] = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            # Skip comments and headers
            if line.startswith('#') or line.startswith('rsid'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                rsid = parts[0]
                allele1 = parts[3]
                allele2 = parts[4]
                
                # Skip invalid/missing genotypes
                if allele1 in ('0', '-', '') or allele2 in ('0', '-', ''):
                    continue
                
                genotypes[rsid] = f"{allele1}{allele2}"
    
    logger.info(f"Loaded {len(genotypes)} genotypes")
    return genotypes


# =============================================================================
# DATASET INITIALIZATION
# =============================================================================

def ensure_datasets_downloaded() -> Dict[str, bool]:
    """
    Ensure all datasets are downloaded and initialized.
    
    Returns:
        Dict mapping dataset name to download success status
    """
    logger.info("Checking dataset status...")
    
    status = get_dataset_status()
    results = {}
    
    datasets_to_init = [
        ('1000genomes', ThousandGenomes),
        ('hgdp', HGDP),
        ('sgdp', SGDPDataset),
        ('ancient_dna', AncientDNADataset),
        ('gnomad', GnomAD),
        ('clinvar', ClinVar),
        ('pharmgkb', PharmGKB),
        ('pgs_catalog', PGSCatalog),
        ('gwas_catalog', GWASCatalog),
    ]
    
    for name, cls in datasets_to_init:
        try:
            ds = cls()
            if not ds.is_downloaded:
                logger.info(f"Downloading {name}...")
                success = ds.download()
                results[name] = success
            else:
                results[name] = True
                logger.info(f"{name}: already downloaded")
        except Exception as e:
            logger.error(f"Error initializing {name}: {e}")
            results[name] = False
    
    return results


# =============================================================================
# POPULATION ANALYSIS
# =============================================================================

def analyze_1kg_populations(
    genotypes: Dict[str, str]
) -> Tuple[List[PopulationResult], Dict[str, float]]:
    """
    Compare genotypes to 1000 Genomes populations.
    
    Args:
        genotypes: Dict mapping rsID to genotype
        
    Returns:
        Tuple of (top populations list, superpopulation breakdown)
    """
    logger.info("Analyzing 1000 Genomes population similarity...")
    
    tg = ThousandGenomes()
    if not tg.is_downloaded:
        tg.download()
    
    similarities = tg.calculate_population_similarity(genotypes)
    superpop_scores = tg.get_superpopulation_summary(similarities)
    
    # Get top populations
    sorted_pops = sorted(similarities.items(), key=lambda x: x[1], reverse=True)
    
    top_pops: List[PopulationResult] = []
    for pop, score in sorted_pops[:10]:
        pop_info = THOUSAND_GENOMES_POPULATIONS.get(pop, {})
        top_pops.append({
            "population": pop,
            "name": pop_info.get("name", pop),
            "region": pop_info.get("superpop", "Unknown"),
            "similarity": round(score * 100, 2),
        })
    
    # Normalize superpopulation scores to percentages
    superpop_pct = {
        SUPERPOPULATIONS.get(k, k): round(v * 100, 2)
        for k, v in superpop_scores.items()
    }
    
    return top_pops, superpop_pct


def analyze_hgdp_populations(
    genotypes: Dict[str, str]
) -> Tuple[List[PopulationResult], Dict[str, float]]:
    """
    Compare genotypes to HGDP populations.
    
    Args:
        genotypes: Dict mapping rsID to genotype
        
    Returns:
        Tuple of (top populations list, region breakdown)
    """
    logger.info("Analyzing HGDP population similarity...")
    
    hgdp = HGDP()
    if not hgdp.is_downloaded:
        hgdp.download()
    
    similarities = hgdp.calculate_population_similarity(genotypes)
    region_scores = hgdp.get_region_summary(similarities)
    
    sorted_pops = sorted(similarities.items(), key=lambda x: x[1], reverse=True)
    
    top_pops: List[PopulationResult] = []
    for pop, score in sorted_pops[:10]:
        pop_info = HGDP_POPULATIONS.get(pop, {})
        top_pops.append({
            "population": pop,
            "name": pop,
            "region": pop_info.get("region", "Unknown"),
            "similarity": round(score * 100, 2),
        })
    
    region_pct = {k: round(v * 100, 2) for k, v in region_scores.items()}
    
    return top_pops, region_pct


def analyze_sgdp_populations(
    genotypes: Dict[str, str]
) -> Tuple[List[PopulationResult], Dict[str, float]]:
    """
    Compare genotypes to SGDP populations.
    
    Args:
        genotypes: Dict mapping rsID to genotype
        
    Returns:
        Tuple of (top populations list, region breakdown)
    """
    logger.info("Analyzing SGDP population similarity...")
    
    try:
        sgdp = SGDPDataset()
        if not sgdp.is_downloaded:
            sgdp.download()
        
        similarities = sgdp.calculate_population_similarity(genotypes)
        
        sorted_pops = sorted(similarities.items(), key=lambda x: x[1], reverse=True)
        
        top_pops: List[PopulationResult] = []
        region_totals: Dict[str, float] = {}
        
        # Import the population data
        from datasets.sgdp import SGDP_POPULATIONS, ALL_SGDP_POPULATIONS
        
        for pop_code, score in sorted_pops[:10]:
            # Find region and name for this population
            region = "Unknown"
            name = pop_code
            for reg, pops in SGDP_POPULATIONS.items():
                for p in pops:
                    if p.code == pop_code:
                        region = reg
                        name = p.name
                        break
            
            top_pops.append({
                "population": pop_code,
                "name": name,
                "region": region,
                "similarity": round(score * 100, 2),
            })
        
        # Aggregate by region
        for pop_code, score in similarities.items():
            for reg, pops in SGDP_POPULATIONS.items():
                if any(p.code == pop_code for p in pops):
                    region_totals[reg] = region_totals.get(reg, 0) + score
                    break
        
        region_pct = {k: round(v * 100, 2) for k, v in region_totals.items()}
        
        return top_pops, region_pct
        
    except Exception as e:
        logger.warning(f"SGDP analysis error: {e}")
        import traceback
        traceback.print_exc()
        return [], {}


# =============================================================================
# ANCIENT ANCESTRY ANALYSIS
# =============================================================================

def analyze_ancient_ancestry(
    genotypes: Dict[str, str]
) -> Tuple[List[AncientSignal], int, int]:
    """
    Analyze ancient ancestry signals.
    
    Args:
        genotypes: Dict mapping rsID to genotype
        
    Returns:
        Tuple of (ancient signals list, neanderthal markers found, total neanderthal markers)
    """
    logger.info("Analyzing ancient ancestry signals...")
    
    try:
        ancient = AncientDNADataset()
        if not ancient.is_downloaded:
            ancient.download()
        
        # Use the dataset's analyze method directly
        results = ancient.analyze_ancestral_signals(genotypes)
        
        signals: List[AncientSignal] = []
        neanderthal_found = 0
        neanderthal_total = 0
        
        for population, signal in results.items():
            signals.append({
                "population": signal.population,
                "signal_strength": signal.signal_strength,
                "markers_matched": signal.markers_found,
                "total_markers": signal.markers_tested,
                "traits": signal.phenotypic_evidence,
            })
            
            if 'neanderthal' in population.lower():
                neanderthal_found = signal.markers_found
                neanderthal_total = signal.markers_tested
        
        return signals, neanderthal_found, neanderthal_total
        
    except Exception as e:
        logger.warning(f"Ancient ancestry analysis error: {e}")
        import traceback
        traceback.print_exc()
        return [], 0, 0


# =============================================================================
# CLINICAL ANALYSIS
# =============================================================================

def analyze_clinvar(
    genotypes: Dict[str, str]
) -> Tuple[List[ClinicalVariant], int, int]:
    """
    Check variants against ClinVar for clinical significance.
    
    Args:
        genotypes: Dict mapping rsID to genotype
        
    Returns:
        Tuple of (clinical findings list, pathogenic count, VUS count)
    """
    logger.info("Analyzing ClinVar clinical variants...")
    
    clinvar = ClinVar()
    if not clinvar.is_downloaded:
        clinvar.download()
    
    findings: List[ClinicalVariant] = []
    pathogenic = 0
    vus = 0
    
    clinvar_results = clinvar.check_user_variants(genotypes)
    
    for rsid, annotation in clinvar_results.items():
        finding: ClinicalVariant = {
            "rsid": rsid,
            "gene": annotation.gene,
            "significance": annotation.clinical_significance,
            "conditions": annotation.conditions,
            "genotype": genotypes.get(rsid, ""),
            "pmids": annotation.pmids,
        }
        findings.append(finding)
        
        if annotation.is_pathogenic:
            pathogenic += 1
        elif "uncertain" in annotation.clinical_significance.lower():
            vus += 1
    
    return findings, pathogenic, vus


def analyze_pharmacogenomics(
    genotypes: Dict[str, str]
) -> Tuple[List[PharmacogenomicResult], int]:
    """
    Analyze pharmacogenomic variants.
    
    Args:
        genotypes: Dict mapping rsID to genotype
        
    Returns:
        Tuple of (pharmacogenomic results list, actionable count)
    """
    logger.info("Analyzing pharmacogenomics...")
    
    pgkb = PharmGKB()
    if not pgkb.is_downloaded:
        pgkb.download()
    
    results: List[PharmacogenomicResult] = []
    actionable = 0
    
    # Key pharmacogenes to analyze
    key_genes = ['CYP2D6', 'CYP2C19', 'CYP2C9', 'DPYD', 'TPMT', 'SLCO1B1', 'VKORC1']
    
    for gene in key_genes:
        try:
            activity_score, phenotype = pgkb.calculate_activity_score(gene, genotypes)
            
            # Get drug interactions for this gene
            interactions = pgkb.get_drug_interactions(gene)
            drugs = list(set(i.drug for i in interactions))
            
            # Get specific recommendation if non-normal
            recommendation = "Standard dosing appropriate"
            if "poor" in phenotype.lower() or "ultra" in phenotype.lower():
                actionable += 1
                # Find specific recommendation
                for interaction in interactions:
                    if phenotype.lower() in interaction.phenotype.lower():
                        recommendation = interaction.recommendation
                        break
            
            results.append({
                "gene": gene,
                "phenotype": phenotype,
                "activity_score": round(activity_score, 2),
                "drugs_affected": drugs[:5],  # Top 5 drugs
                "recommendation": recommendation,
            })
        except Exception as e:
            logger.warning(f"Error analyzing {gene}: {e}")
    
    return results, actionable


# =============================================================================
# PRS ANALYSIS
# =============================================================================

def analyze_polygenic_scores(
    genotypes: Dict[str, str]
) -> Tuple[List[PRSResult], List[str]]:
    """
    Calculate polygenic risk scores for available conditions.
    
    Args:
        genotypes: Dict mapping rsID to genotype
        
    Returns:
        Tuple of (PRS results list, high-risk conditions list)
    """
    logger.info("Calculating polygenic risk scores...")
    
    pgs = PGSCatalog()
    if not pgs.is_downloaded:
        pgs.download()
    
    results: List[PRSResult] = []
    high_risk: List[str] = []
    
    all_prs = pgs.calculate_all_prs(genotypes)
    
    for trait, prs_result in all_prs.items():
        result: PRSResult = {
            "trait": prs_result.trait,
            "raw_score": round(prs_result.raw_score, 3),
            "percentile": round(prs_result.percentile, 1),
            "risk_category": prs_result.risk_category,
            "interpretation": prs_result.interpretation,
            "pmid": prs_result.pmid,
        }
        results.append(result)
        
        if prs_result.risk_category in ["High", "Elevated"]:
            high_risk.append(trait)
    
    return results, high_risk


# =============================================================================
# GWAS ANALYSIS
# =============================================================================

def analyze_gwas_traits(
    genotypes: Dict[str, str]
) -> Tuple[List[GWASFinding], List[str]]:
    """
    Check variants against GWAS Catalog for trait associations.
    
    Args:
        genotypes: Dict mapping rsID to genotype
        
    Returns:
        Tuple of (GWAS findings list, notable traits list)
    """
    logger.info("Analyzing GWAS trait associations...")
    
    gwas = GWASCatalog()
    if not gwas.is_downloaded:
        gwas.download()
    
    findings: List[GWASFinding] = []
    notable_traits: List[str] = []
    
    user_variants = gwas.check_user_variants(genotypes)
    
    for rsid, associations in user_variants.items():
        user_geno = genotypes.get(rsid, "")
        
        for assoc in associations:
            # Check if user carries risk allele
            risk_allele = assoc.risk_allele.upper()
            has_risk = risk_allele in user_geno.upper()
            
            finding: GWASFinding = {
                "rsid": rsid,
                "trait": assoc.trait,
                "risk_allele": assoc.risk_allele,
                "user_genotype": user_geno,
                "effect_direction": "carries risk allele" if has_risk else "does not carry risk allele",
                "odds_ratio": assoc.odds_ratio,
                "pmid": assoc.pmid,
            }
            findings.append(finding)
            
            # Track notable findings (high OR + carries risk)
            if has_risk and assoc.odds_ratio and assoc.odds_ratio > 1.5:
                notable_traits.append(assoc.trait)
    
    return findings, list(set(notable_traits))


# =============================================================================
# MAIN ANALYSIS FUNCTION
# =============================================================================

def run_comprehensive_analysis(
    dna_filepath: str,
    output_dir: Optional[str] = None
) -> ComprehensiveAnalysisResult:
    """
    Run complete analysis on DNA file using all 9 datasets.
    
    Args:
        dna_filepath: Path to the DNA data file
        output_dir: Optional output directory for reports
        
    Returns:
        ComprehensiveAnalysisResult with all findings
    """
    logger.info("=" * 60)
    logger.info("STARTING COMPREHENSIVE DNA ANALYSIS")
    logger.info("=" * 60)
    
    start_time = datetime.now()
    
    # Initialize result
    result = ComprehensiveAnalysisResult(
        analysis_date=start_time.isoformat(),
        dna_file=dna_filepath,
    )
    
    # Step 1: Ensure all datasets are ready
    logger.info("\n[1/8] Initializing datasets...")
    dataset_status = ensure_datasets_downloaded()
    logger.info(f"Dataset status: {dataset_status}")
    
    # Step 2: Load DNA data
    logger.info("\n[2/8] Loading DNA data...")
    genotypes = load_dna_file(dna_filepath)
    result.total_snps = len(genotypes)
    
    # Step 3: 1000 Genomes population analysis
    logger.info("\n[3/8] Analyzing 1000 Genomes populations...")
    top_1kg, superpop = analyze_1kg_populations(genotypes)
    result.top_1kg_populations = top_1kg
    result.superpopulation_breakdown = superpop
    
    # Step 4: HGDP analysis
    logger.info("\n[4/8] Analyzing HGDP populations...")
    top_hgdp, hgdp_regions = analyze_hgdp_populations(genotypes)
    result.top_hgdp_populations = top_hgdp
    result.hgdp_region_breakdown = hgdp_regions
    
    # Step 5: SGDP analysis
    logger.info("\n[5/8] Analyzing SGDP populations...")
    top_sgdp, sgdp_regions = analyze_sgdp_populations(genotypes)
    result.top_sgdp_populations = top_sgdp
    result.sgdp_region_breakdown = sgdp_regions
    
    # Step 6: Ancient ancestry
    logger.info("\n[6/8] Analyzing ancient ancestry...")
    ancient, neand_found, neand_total = analyze_ancient_ancestry(genotypes)
    result.ancient_signals = ancient
    result.neanderthal_markers = neand_found
    result.neanderthal_total = neand_total
    
    # Step 7: Clinical analysis
    logger.info("\n[7/8] Analyzing clinical variants...")
    clinvar_findings, pathogenic, vus = analyze_clinvar(genotypes)
    result.clinvar_findings = clinvar_findings
    result.pathogenic_count = pathogenic
    result.vus_count = vus
    
    pgx_results, actionable_pgx = analyze_pharmacogenomics(genotypes)
    result.pharmacogenomics = pgx_results
    result.actionable_pgx = actionable_pgx
    
    # Step 8: Risk assessment
    logger.info("\n[8/8] Calculating risk scores and trait associations...")
    prs_results, high_risk = analyze_polygenic_scores(genotypes)
    result.prs_results = prs_results
    result.high_risk_conditions = high_risk
    
    gwas_findings, notable_traits = analyze_gwas_traits(genotypes)
    result.gwas_findings = gwas_findings
    result.notable_traits = notable_traits
    
    # Calculate analysis stats
    result.snps_analyzed = len(set(
        [f["rsid"] for f in clinvar_findings] +
        [f["rsid"] for f in gwas_findings]
    ))
    
    # Timing
    elapsed = (datetime.now() - start_time).total_seconds()
    logger.info(f"\nAnalysis completed in {elapsed:.1f} seconds")
    
    # Save results
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save JSON
        json_path = output_path / "comprehensive_analysis.json"
        with open(json_path, 'w') as f:
            json.dump(asdict(result), f, indent=2, default=str)
        logger.info(f"Saved JSON results to {json_path}")
    
    return result


# =============================================================================
# REPORT GENERATION
# =============================================================================

def generate_text_report(result: ComprehensiveAnalysisResult) -> str:
    """Generate human-readable text report."""
    
    lines = [
        "=" * 70,
        "COMPREHENSIVE DNA ANALYSIS REPORT",
        "=" * 70,
        "",
        f"Analysis Date: {result.analysis_date}",
        f"DNA File: {result.dna_file}",
        f"Total SNPs: {result.total_snps:,}",
        "",
        "-" * 70,
        "POPULATION ANCESTRY",
        "-" * 70,
        "",
        "1000 Genomes - Superpopulation Breakdown:",
    ]
    
    for superpop, pct in sorted(result.superpopulation_breakdown.items(), 
                                 key=lambda x: x[1], reverse=True):
        lines.append(f"  • {superpop}: {pct:.1f}%")
    
    lines.extend([
        "",
        "Top Similar Populations (1000 Genomes):",
    ])
    
    for pop in result.top_1kg_populations[:5]:
        lines.append(f"  • {pop['name']} ({pop['population']}): {pop['similarity']:.1f}%")
    
    lines.extend([
        "",
        "Top Similar Populations (HGDP):",
    ])
    
    for pop in result.top_hgdp_populations[:5]:
        lines.append(f"  • {pop['name']}: {pop['similarity']:.1f}%")
    
    lines.extend([
        "",
        "-" * 70,
        "ANCIENT ANCESTRY SIGNALS",
        "-" * 70,
        "",
    ])
    
    if result.ancient_signals:
        for signal in result.ancient_signals:
            lines.append(f"• {signal['population']}: {signal['signal_strength'].upper()}")
            lines.append(f"  Markers: {signal['markers_matched']}/{signal['total_markers']}")
            if signal['traits']:
                lines.append(f"  Associated traits: {', '.join(signal['traits'][:3])}")
    else:
        lines.append("  (No ancient signals analyzed)")
    
    lines.extend([
        "",
        f"Neanderthal Markers: {result.neanderthal_markers}/{result.neanderthal_total}",
        "",
        "-" * 70,
        "CLINICAL FINDINGS",
        "-" * 70,
        "",
        f"Pathogenic variants: {result.pathogenic_count}",
        f"Variants of uncertain significance: {result.vus_count}",
        "",
    ])
    
    if result.clinvar_findings:
        for finding in result.clinvar_findings[:5]:
            lines.append(f"• {finding['rsid']} ({finding['gene']}): {finding['significance']}")
            if finding['conditions']:
                lines.append(f"  Conditions: {', '.join(finding['conditions'][:2])}")
    
    lines.extend([
        "",
        "-" * 70,
        "PHARMACOGENOMICS",
        "-" * 70,
        "",
        f"Actionable findings: {result.actionable_pgx}",
        "",
    ])
    
    for pgx in result.pharmacogenomics:
        lines.append(f"• {pgx['gene']}: {pgx['phenotype']} (activity: {pgx['activity_score']})")
        if pgx['drugs_affected']:
            lines.append(f"  Affects: {', '.join(pgx['drugs_affected'][:3])}")
    
    lines.extend([
        "",
        "-" * 70,
        "POLYGENIC RISK SCORES",
        "-" * 70,
        "",
    ])
    
    for prs in result.prs_results:
        status = "⚠️" if prs['risk_category'] in ['High', 'Elevated'] else "✓"
        lines.append(f"{status} {prs['trait']}: {prs['risk_category']} (percentile: {prs['percentile']:.0f})")
    
    if result.high_risk_conditions:
        lines.extend([
            "",
            "⚠️ Elevated Risk Conditions:",
        ])
        for condition in result.high_risk_conditions:
            lines.append(f"  • {condition}")
    
    lines.extend([
        "",
        "-" * 70,
        "TRAIT ASSOCIATIONS (GWAS)",
        "-" * 70,
        "",
    ])
    
    if result.notable_traits:
        lines.append("Notable trait associations:")
        for trait in result.notable_traits[:10]:
            lines.append(f"  • {trait}")
    
    lines.extend([
        "",
        "=" * 70,
        "END OF REPORT",
        "=" * 70,
    ])
    
    return "\n".join(lines)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    # Default paths
    dna_file = Path.home() / "dna-analysis" / "raw_data.txt"
    output_dir = Path.home() / "dna-analysis" / "reports" / "overnight_2026-02-07"
    
    # Run analysis
    result = run_comprehensive_analysis(str(dna_file), str(output_dir))
    
    # Generate text report
    report = generate_text_report(result)
    report_path = output_dir / "comprehensive_report.txt"
    with open(report_path, 'w') as f:
        f.write(report)
    
    logger.info(f"\nReport saved to: {report_path}")
    logger.info("\n" + report)
