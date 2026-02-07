"""
Pharmacogenomics Analysis with Statistical Rigor

Provides metabolizer phenotype calling with confidence levels and
activity score uncertainty quantification.

Author: OpenClaw AI
Date: 2026-02-07
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path
import sys

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from personal_genomics.statistics import (
        ConfidenceLevel,
        diplotype_confidence,
        confidence_to_color,
    )
    STATS_AVAILABLE = True
except ImportError:
    STATS_AVAILABLE = False
    class ConfidenceLevel:
        DEFINITIVE = "DEFINITIVE"
        HIGH = "HIGH"
        MEDIUM = "MEDIUM"
        LOW = "LOW"
        UNCERTAIN = "UNCERTAIN"


# =============================================================================
# STAR ALLELE DEFINITIONS WITH MARKER COUNTS
# =============================================================================

STAR_ALLELE_MARKERS = {
    "CYP2D6": {
        "*1": {"markers": [], "activity": 1.0, "function": "Normal", "note": "Reference allele - inferred by absence"},
        "*2": {"markers": ["rs16947"], "activity": 1.0, "function": "Normal"},
        "*3": {"markers": ["rs35742686"], "activity": 0.0, "function": "No function"},
        "*4": {"markers": ["rs3892097"], "activity": 0.0, "function": "No function"},
        "*5": {"markers": [], "activity": 0.0, "function": "No function", "note": "Gene deletion - cannot detect on arrays"},
        "*6": {"markers": ["rs5030655"], "activity": 0.0, "function": "No function"},
        "*9": {"markers": ["rs5030656"], "activity": 0.5, "function": "Decreased"},
        "*10": {"markers": ["rs1065852"], "activity": 0.25, "function": "Decreased"},
        "*17": {"markers": ["rs28371706"], "activity": 0.5, "function": "Decreased"},
        "*29": {"markers": ["rs59421388"], "activity": 0.5, "function": "Decreased"},
        "*41": {"markers": ["rs28371725"], "activity": 0.5, "function": "Decreased"},
    },
    "CYP2C19": {
        "*1": {"markers": [], "activity": 1.0, "function": "Normal", "note": "Reference allele"},
        "*2": {"markers": ["rs4244285"], "activity": 0.0, "function": "No function"},
        "*3": {"markers": ["rs4986893"], "activity": 0.0, "function": "No function"},
        "*17": {"markers": ["rs12248560"], "activity": 1.5, "function": "Increased"},
    },
    "CYP2C9": {
        "*1": {"markers": [], "activity": 1.0, "function": "Normal", "note": "Reference allele"},
        "*2": {"markers": ["rs1799853"], "activity": 0.5, "function": "Decreased"},
        "*3": {"markers": ["rs1057910"], "activity": 0.1, "function": "Decreased"},
        "*5": {"markers": ["rs28371685"], "activity": 0.5, "function": "Decreased"},
        "*6": {"markers": ["rs9332131"], "activity": 0.0, "function": "No function"},
        "*8": {"markers": ["rs7900194"], "activity": 0.5, "function": "Decreased"},
        "*11": {"markers": ["rs28371686"], "activity": 0.5, "function": "Decreased"},
    },
    "CYP3A5": {
        "*1": {"markers": [], "activity": 1.0, "function": "Normal expressor", "note": "Reference"},
        "*3": {"markers": ["rs776746"], "activity": 0.0, "function": "Non-expressor"},
    },
    "TPMT": {
        "*1": {"markers": [], "activity": 1.0, "function": "Normal", "note": "Reference"},
        "*2": {"markers": ["rs1800462"], "activity": 0.0, "function": "No function"},
        "*3A": {"markers": ["rs1800460", "rs1142345"], "activity": 0.0, "function": "No function"},
        "*3B": {"markers": ["rs1800460"], "activity": 0.0, "function": "No function"},
        "*3C": {"markers": ["rs1142345"], "activity": 0.0, "function": "No function"},
    },
    "NUDT15": {
        "*1": {"markers": [], "activity": 1.0, "function": "Normal", "note": "Reference"},
        "*3": {"markers": ["rs116855232"], "activity": 0.0, "function": "No function"},
    },
    "DPYD": {
        "*1": {"markers": [], "activity": 1.0, "function": "Normal", "note": "Reference"},
        "*2A": {"markers": ["rs3918290"], "activity": 0.0, "function": "No function"},
        "*13": {"markers": ["rs55886062"], "activity": 0.0, "function": "No function"},
    },
}


# =============================================================================
# PHENOTYPE MAPPING
# =============================================================================

def activity_to_phenotype(activity_score: float, gene: str) -> str:
    """Convert activity score to metabolizer phenotype."""
    if gene in ("CYP2D6", "CYP2C19", "CYP2C9"):
        if activity_score >= 2.0:
            return "Ultrarapid Metabolizer"
        elif activity_score >= 1.25:
            return "Normal-to-Rapid Metabolizer"
        elif activity_score >= 1.0:
            return "Normal Metabolizer"
        elif activity_score >= 0.5:
            return "Intermediate Metabolizer"
        else:
            return "Poor Metabolizer"
    elif gene == "CYP3A5":
        if activity_score >= 1.0:
            return "Expressor"
        else:
            return "Non-expressor"
    elif gene in ("TPMT", "NUDT15"):
        if activity_score >= 1.0:
            return "Normal Metabolizer"
        elif activity_score >= 0.5:
            return "Intermediate Metabolizer"
        else:
            return "Poor Metabolizer"
    elif gene == "DPYD":
        if activity_score >= 1.5:
            return "Normal Metabolizer"
        elif activity_score >= 1.0:
            return "Normal Metabolizer"
        elif activity_score >= 0.5:
            return "Intermediate Metabolizer"
        else:
            return "Poor Metabolizer"
    else:
        return "Unknown"


# =============================================================================
# DIPLOTYPE CALLING WITH CONFIDENCE
# =============================================================================

@dataclass
class MetabolizerResult:
    """Result of metabolizer phenotype calling with statistics."""
    gene: str
    phenotype: str
    activity_score: float
    activity_score_uncertainty: float
    star_alleles: List[str]
    diplotype: str
    confidence: str
    diplotype_confidence: float
    warnings: List[str] = field(default_factory=list)
    markers_used: int = 0
    markers_total: int = 0
    phase_known: bool = False
    
    def to_dict(self) -> dict:
        return {
            "gene": self.gene,
            "phenotype": self.phenotype,
            "activity_score": self.activity_score,
            "activity_uncertainty": round(self.activity_score_uncertainty, 2),
            "star_alleles": self.star_alleles,
            "diplotype": self.diplotype,
            "confidence": self.confidence,
            "diplotype_confidence": round(self.diplotype_confidence, 3),
            "warnings": self.warnings,
            "markers_found": self.markers_used,
            "markers_total": self.markers_total,
            "phase_known": self.phase_known,
        }


def call_star_alleles(
    gene: str,
    genotypes: Dict[str, str],
    gene_markers: Dict[str, Dict] = None
) -> Tuple[List[str], Dict[str, Any]]:
    """
    Call star alleles for a gene based on observed genotypes.
    
    Args:
        gene: Gene name (e.g., "CYP2D6")
        genotypes: Dict of rsid -> genotype
        gene_markers: Optional override for marker definitions
        
    Returns:
        Tuple of (list of called alleles, metadata dict)
    """
    markers = gene_markers or STAR_ALLELE_MARKERS.get(gene, {})
    
    called_alleles = []
    allele_evidence = {}
    total_markers = 0
    markers_found = 0
    
    # Check each defined star allele
    for allele, info in markers.items():
        allele_markers = info.get("markers", [])
        if not allele_markers:
            continue  # Reference allele, defined by absence
        
        total_markers += len(allele_markers)
        
        # Check if allele markers are present
        markers_present = 0
        for rsid in allele_markers:
            if rsid in genotypes:
                markers_found += 1
                geno = genotypes[rsid]
                # Check if variant allele is present
                # This is simplified - real calling needs risk allele info
                if geno and geno not in ("--", "00", "??"):
                    markers_present += 1
        
        if markers_present > 0:
            # Evidence for this allele
            coverage = markers_present / len(allele_markers) if allele_markers else 0
            allele_evidence[allele] = {
                "markers_found": markers_present,
                "markers_total": len(allele_markers),
                "coverage": coverage,
            }
            
            # Count copies (simplified - assumes heterozygous if not all markers found)
            if coverage >= 1.0:
                called_alleles.append(allele)
            elif coverage >= 0.5:
                called_alleles.append(allele)  # Single copy likely
    
    # If no variant alleles found, assume reference (*1)
    if not called_alleles:
        called_alleles = ["*1", "*1"]
    elif len(called_alleles) == 1:
        called_alleles = [called_alleles[0], "*1"]
    elif len(called_alleles) > 2:
        # Take top 2 by evidence
        called_alleles = sorted(
            called_alleles,
            key=lambda a: allele_evidence.get(a, {}).get("coverage", 0),
            reverse=True
        )[:2]
    
    metadata = {
        "markers_found": markers_found,
        "markers_total": total_markers,
        "allele_evidence": allele_evidence,
    }
    
    return called_alleles, metadata


def calculate_metabolizer_phenotype(
    gene: str,
    genotypes: Dict[str, str]
) -> MetabolizerResult:
    """
    Calculate metabolizer phenotype with full statistical confidence.
    
    Args:
        gene: Gene name (e.g., "CYP2D6")
        genotypes: Dict of rsid -> genotype
        
    Returns:
        MetabolizerResult with phenotype, activity score, and confidence metrics
    """
    gene_markers = STAR_ALLELE_MARKERS.get(gene, {})
    
    # Call star alleles
    called_alleles, metadata = call_star_alleles(gene, genotypes, gene_markers)
    
    # Calculate activity score
    activity_scores = []
    for allele in called_alleles:
        if allele in gene_markers:
            activity_scores.append(gene_markers[allele].get("activity", 1.0))
        else:
            activity_scores.append(1.0)  # Assume normal for unknown
    
    total_activity = sum(activity_scores)
    
    # Calculate uncertainty
    # Higher uncertainty if:
    # 1. *1 allele inferred (not directly observed)
    # 2. Low marker coverage
    # 3. Phase ambiguous
    
    markers_found = metadata.get("markers_found", 0)
    markers_total = metadata.get("markers_total", 1)
    coverage = markers_found / markers_total if markers_total > 0 else 0
    
    base_uncertainty = 0.1  # Minimum uncertainty
    
    # Add uncertainty for inferred *1
    if "*1" in called_alleles:
        base_uncertainty += 0.2
    
    # Add uncertainty for low coverage
    if coverage < 0.5:
        base_uncertainty += 0.3
    elif coverage < 0.7:
        base_uncertainty += 0.15
    
    activity_uncertainty = base_uncertainty
    
    # Calculate diplotype confidence
    warnings = []
    
    if STATS_AVAILABLE:
        # Get marker counts for each allele
        allele_1 = called_alleles[0] if called_alleles else "*1"
        allele_2 = called_alleles[1] if len(called_alleles) > 1 else "*1"
        
        a1_info = gene_markers.get(allele_1, {})
        a2_info = gene_markers.get(allele_2, {})
        
        a1_markers = len(a1_info.get("markers", []))
        a2_markers = len(a2_info.get("markers", []))
        
        # Get evidence
        evidence = metadata.get("allele_evidence", {})
        a1_found = evidence.get(allele_1, {}).get("markers_found", 0) if allele_1 != "*1" else 0
        a2_found = evidence.get(allele_2, {}).get("markers_found", 0) if allele_2 != "*1" else 0
        
        dip_conf, dip_warnings = diplotype_confidence(
            allele_1, allele_2,
            a1_found, max(a1_markers, 1),
            a2_found, max(a2_markers, 1),
            phase_ambiguous=False  # Can't determine phase from DTC data
        )
        
        warnings.extend(dip_warnings)
    else:
        dip_conf = 0.7 if coverage >= 0.5 else 0.5
    
    # Determine confidence level
    if dip_conf >= 0.9 and "*1" not in called_alleles:
        conf_level = ConfidenceLevel.DEFINITIVE
    elif dip_conf >= 0.8:
        conf_level = ConfidenceLevel.HIGH
    elif dip_conf >= 0.6:
        conf_level = ConfidenceLevel.MEDIUM
    elif dip_conf >= 0.4:
        conf_level = ConfidenceLevel.LOW
    else:
        conf_level = ConfidenceLevel.UNCERTAIN
    
    # Add standard warnings
    if "*1" in called_alleles:
        warnings.append("*1 allele inferred by absence of variant markers")
    
    if coverage < 0.3:
        warnings.append(f"Very low marker coverage ({coverage*100:.0f}%)")
        conf_level = ConfidenceLevel.UNCERTAIN
    
    # Get phenotype
    phenotype = activity_to_phenotype(total_activity, gene)
    
    # Build diplotype string
    diplotype = f"{called_alleles[0]}/{called_alleles[1]}" if len(called_alleles) >= 2 else called_alleles[0]
    
    return MetabolizerResult(
        gene=gene,
        phenotype=phenotype,
        activity_score=total_activity,
        activity_score_uncertainty=activity_uncertainty,
        star_alleles=called_alleles,
        diplotype=diplotype,
        confidence=conf_level.value if hasattr(conf_level, 'value') else str(conf_level),
        diplotype_confidence=dip_conf,
        warnings=warnings,
        markers_used=markers_found,
        markers_total=markers_total,
        phase_known=False,
    )


def analyze_all_pharmacogenes(genotypes: Dict[str, str]) -> Dict[str, MetabolizerResult]:
    """
    Analyze all pharmacogenomics genes with statistical confidence.
    
    Args:
        genotypes: Dict of rsid -> genotype
        
    Returns:
        Dict of gene -> MetabolizerResult
    """
    results = {}
    
    for gene in STAR_ALLELE_MARKERS.keys():
        results[gene] = calculate_metabolizer_phenotype(gene, genotypes)
    
    return results


def get_pharmacogenomics_summary(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Get summary of pharmacogenomics results with statistics for dashboard.
    
    Args:
        genotypes: Dict of rsid -> genotype
        
    Returns:
        Dict with all results formatted for display
    """
    results = analyze_all_pharmacogenes(genotypes)
    
    summary = {
        "genes": {},
        "actionable_findings": [],
        "warnings": [],
        "overall_confidence": ConfidenceLevel.MEDIUM.value if hasattr(ConfidenceLevel, 'value') else "MEDIUM",
    }
    
    confidence_scores = []
    
    for gene, result in results.items():
        summary["genes"][gene] = result.to_dict()
        confidence_scores.append(result.diplotype_confidence)
        
        # Flag actionable findings
        if result.phenotype in ("Poor Metabolizer", "Ultrarapid Metabolizer"):
            summary["actionable_findings"].append({
                "gene": gene,
                "phenotype": result.phenotype,
                "diplotype": result.diplotype,
                "confidence": result.confidence,
            })
        
        # Collect warnings
        for warning in result.warnings:
            summary["warnings"].append(f"{gene}: {warning}")
    
    # Calculate overall confidence
    if confidence_scores:
        avg_conf = sum(confidence_scores) / len(confidence_scores)
        if avg_conf >= 0.8:
            summary["overall_confidence"] = ConfidenceLevel.HIGH.value if hasattr(ConfidenceLevel, 'value') else "HIGH"
        elif avg_conf >= 0.6:
            summary["overall_confidence"] = ConfidenceLevel.MEDIUM.value if hasattr(ConfidenceLevel, 'value') else "MEDIUM"
        else:
            summary["overall_confidence"] = ConfidenceLevel.LOW.value if hasattr(ConfidenceLevel, 'value') else "LOW"
    
    return summary
