"""
Data Quality Assessment Module for Personal Genomics

Provides comprehensive quality metrics for genotyping data including:
- Overall chip quality scoring
- Call rate analysis by chromosome
- Marker coverage by category
- Missing critical markers flagging
- Platform detection confidence

Author: OpenClaw AI
Date: 2026-02-07
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Set, Tuple, Any
from collections import defaultdict

from .statistics import ConfidenceLevel, StatisticalResult, wilson_score_interval


# =============================================================================
# TYPE DEFINITIONS
# =============================================================================

class QualityGrade(Enum):
    """Overall quality grade for genotyping data."""
    EXCELLENT = "A"  # >95% call rate, all critical markers present
    GOOD = "B"       # >90% call rate, most critical markers present
    FAIR = "C"       # >80% call rate, some critical markers missing
    POOR = "D"       # <80% call rate or many critical markers missing
    FAIL = "F"       # Data unsuitable for analysis


@dataclass
class ChromosomeQuality:
    """Quality metrics for a single chromosome."""
    chromosome: str
    total_snps: int
    called_snps: int
    no_calls: int
    call_rate: float
    expected_min: int
    expected_max: int
    coverage_status: str  # "normal", "low", "high", "missing"
    
    def to_dict(self) -> dict:
        return {
            "chromosome": self.chromosome,
            "total_snps": self.total_snps,
            "called_snps": self.called_snps,
            "no_calls": self.no_calls,
            "call_rate": round(self.call_rate, 4),
            "expected_range": [self.expected_min, self.expected_max],
            "coverage_status": self.coverage_status
        }


@dataclass
class CategoryCoverage:
    """Coverage metrics for an analysis category."""
    category: str
    markers_in_reference: int
    markers_found: int
    markers_called: int  # Found and have valid genotype
    coverage_rate: float
    critical_missing: List[str] = field(default_factory=list)
    confidence_level: ConfidenceLevel = ConfidenceLevel.UNCERTAIN
    
    def to_dict(self) -> dict:
        return {
            "category": self.category,
            "markers_in_reference": self.markers_in_reference,
            "markers_found": self.markers_found,
            "markers_called": self.markers_called,
            "coverage_rate": round(self.coverage_rate, 4),
            "critical_missing": self.critical_missing,
            "confidence": self.confidence_level.value
        }


@dataclass
class PlatformDetection:
    """Platform/chip detection results."""
    detected_platform: str
    confidence: float
    total_snps: int
    platform_details: Optional[str] = None
    alternative_matches: List[str] = field(default_factory=list)
    signature_snps_found: int = 0
    signature_snps_total: int = 0
    
    def to_dict(self) -> dict:
        return {
            "platform": self.detected_platform,
            "confidence": round(self.confidence, 3),
            "total_snps": self.total_snps,
            "details": self.platform_details,
            "alternatives": self.alternative_matches,
            "signature_snps": f"{self.signature_snps_found}/{self.signature_snps_total}"
        }


@dataclass
class QualityReport:
    """Comprehensive quality assessment report."""
    overall_grade: QualityGrade
    overall_score: float  # 0-100
    total_snps: int
    called_snps: int
    overall_call_rate: float
    chromosome_quality: Dict[str, ChromosomeQuality]
    category_coverage: Dict[str, CategoryCoverage]
    platform: PlatformDetection
    critical_markers_status: Dict[str, str]  # rsid -> "present"|"missing"|"no_call"
    warnings: List[str] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)
    confidence_level: ConfidenceLevel = ConfidenceLevel.MEDIUM
    
    def to_dict(self) -> dict:
        return {
            "overall_grade": self.overall_grade.value,
            "overall_score": round(self.overall_score, 1),
            "total_snps": self.total_snps,
            "called_snps": self.called_snps,
            "overall_call_rate": round(self.overall_call_rate, 4),
            "chromosome_quality": {
                k: v.to_dict() for k, v in self.chromosome_quality.items()
            },
            "category_coverage": {
                k: v.to_dict() for k, v in self.category_coverage.items()
            },
            "platform": self.platform.to_dict(),
            "critical_markers": self.critical_markers_status,
            "warnings": self.warnings,
            "recommendations": self.recommendations,
            "confidence": self.confidence_level.value
        }


# =============================================================================
# CRITICAL MARKERS - Must flag if missing
# =============================================================================

CRITICAL_MARKERS = {
    # Pharmacogenomics - Life-threatening drug interactions
    "rs3918290": {"gene": "DPYD", "reason": "5-FU/capecitabine fatal toxicity"},
    "rs2395029": {"gene": "HLA-B*5701", "reason": "Abacavir hypersensitivity"},
    "rs3909184": {"gene": "HLA-B*1502", "reason": "Carbamazepine SJS/TEN"},
    "rs1142345": {"gene": "TPMT", "reason": "Thiopurine myelosuppression"},
    "rs116855232": {"gene": "NUDT15", "reason": "Thiopurine toxicity (Asian)"},
    
    # Warfarin dosing
    "rs9923231": {"gene": "VKORC1", "reason": "Warfarin dose requirement"},
    "rs1799853": {"gene": "CYP2C9*2", "reason": "Warfarin metabolism"},
    "rs1057910": {"gene": "CYP2C9*3", "reason": "Warfarin metabolism"},
    
    # Common high-impact PGx
    "rs4244285": {"gene": "CYP2C19*2", "reason": "Clopidogrel efficacy"},
    "rs3892097": {"gene": "CYP2D6*4", "reason": "Codeine/tamoxifen metabolism"},
    "rs4149056": {"gene": "SLCO1B1", "reason": "Statin myopathy"},
    
    # Disease risk
    "rs429358": {"gene": "APOE", "reason": "Alzheimer's risk (ε4)"},
    "rs7412": {"gene": "APOE", "reason": "Alzheimer's risk (ε2/ε4)"},
    "rs6025": {"gene": "F5 Leiden", "reason": "VTE/contraceptive risk"},
    "rs1799963": {"gene": "F2", "reason": "Prothrombin mutation"},
    
    # Carrier status
    "rs80357906": {"gene": "BRCA1", "reason": "Hereditary breast cancer"},
    "rs80359550": {"gene": "BRCA2", "reason": "Hereditary breast cancer"},
}


# =============================================================================
# EXPECTED CHROMOSOME COVERAGE
# =============================================================================

CHROMOSOME_EXPECTED = {
    "1": (45000, 75000),
    "2": (45000, 70000),
    "3": (35000, 55000),
    "4": (30000, 50000),
    "5": (30000, 50000),
    "6": (35000, 55000),
    "7": (30000, 45000),
    "8": (25000, 40000),
    "9": (25000, 40000),
    "10": (25000, 40000),
    "11": (25000, 40000),
    "12": (25000, 40000),
    "13": (18000, 30000),
    "14": (17000, 28000),
    "15": (16000, 26000),
    "16": (17000, 28000),
    "17": (16000, 26000),
    "18": (15000, 24000),
    "19": (12000, 22000),
    "20": (13000, 22000),
    "21": (7000, 14000),
    "22": (8000, 14000),
    "X": (20000, 40000),
    "Y": (500, 5000),
    "MT": (100, 1000),
}


# =============================================================================
# PLATFORM SIGNATURES
# =============================================================================

PLATFORM_SIGNATURES = {
    "23andme_v5": {
        "marker_count_range": (630000, 680000),
        "signature_snps": ["rs548049170", "rs9461509", "rs73211851", "rs796052984"],
        "description": "23andMe v5 chip (2017+)"
    },
    "23andme_v4": {
        "marker_count_range": (570000, 620000),
        "signature_snps": ["rs12203592", "rs6548616", "rs4988235"],
        "description": "23andMe v4 chip (2013-2017)"
    },
    "23andme_v3": {
        "marker_count_range": (950000, 1050000),
        "signature_snps": ["rs2032658", "rs12785878"],
        "description": "23andMe v3 chip (2010-2013)"
    },
    "ancestrydna_v2": {
        "marker_count_range": (680000, 750000),
        "signature_snps": ["rs1667394", "rs4778138"],
        "description": "AncestryDNA v2 chip"
    },
    "myheritage": {
        "marker_count_range": (680000, 750000),
        "signature_snps": [],
        "description": "MyHeritage DNA"
    },
    "ftdna": {
        "marker_count_range": (680000, 750000),
        "signature_snps": [],
        "description": "FamilyTreeDNA"
    },
}


# =============================================================================
# CATEGORY MARKER COUNTS (Reference)
# =============================================================================

CATEGORY_MARKER_COUNTS = {
    "pharmacogenomics": {"total": 159, "critical": 15},
    "polygenic_risk_scores": {"total": 277, "critical": 20},
    "carrier_status": {"total": 181, "critical": 10},
    "health_risks": {"total": 233, "critical": 15},
    "traits": {"total": 163, "critical": 5},
    "haplogroups": {"total": 44, "critical": 10},
    "ancestry": {"total": 124, "critical": 0},
    "cancer_panel": {"total": 41, "critical": 10},
    "autoimmune": {"total": 31, "critical": 5},
    "ancient_dna": {"total": 50, "critical": 0},
}


# =============================================================================
# QUALITY ASSESSMENT FUNCTIONS
# =============================================================================

def assess_chromosome_quality(
    genotypes: Dict[str, str],
    chromosome_map: Dict[str, str]  # rsid -> chromosome
) -> Dict[str, ChromosomeQuality]:
    """
    Assess quality metrics per chromosome.
    
    Args:
        genotypes: Dict of rsid -> genotype
        chromosome_map: Dict of rsid -> chromosome number
        
    Returns:
        Dict of chromosome -> ChromosomeQuality
    """
    # Count SNPs per chromosome
    chrom_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: {"total": 0, "called": 0, "nocall": 0})
    
    for rsid, geno in genotypes.items():
        chrom = chromosome_map.get(rsid, "unknown")
        if chrom == "unknown":
            continue
        
        # Normalize chromosome name
        chrom = str(chrom).upper().replace("CHR", "")
        
        chrom_counts[chrom]["total"] += 1
        
        if geno and geno not in ("--", "00", "??", "NC", ""):
            chrom_counts[chrom]["called"] += 1
        else:
            chrom_counts[chrom]["nocall"] += 1
    
    # Build quality reports
    results = {}
    for chrom, counts in chrom_counts.items():
        expected = CHROMOSOME_EXPECTED.get(chrom, (1000, 100000))
        total = counts["total"]
        called = counts["called"]
        nocall = counts["nocall"]
        
        call_rate = called / total if total > 0 else 0
        
        # Determine coverage status
        if total < expected[0] * 0.5:
            status = "low"
        elif total > expected[1] * 1.5:
            status = "high"  # Unusual, might be WGS
        else:
            status = "normal"
        
        results[chrom] = ChromosomeQuality(
            chromosome=chrom,
            total_snps=total,
            called_snps=called,
            no_calls=nocall,
            call_rate=call_rate,
            expected_min=expected[0],
            expected_max=expected[1],
            coverage_status=status
        )
    
    return results


def assess_category_coverage(
    genotypes: Dict[str, str],
    category_markers: Dict[str, Set[str]],  # category -> set of rsids
    critical_per_category: Dict[str, Set[str]] = None
) -> Dict[str, CategoryCoverage]:
    """
    Assess marker coverage for each analysis category.
    
    Args:
        genotypes: Dict of rsid -> genotype
        category_markers: Dict of category -> set of marker rsids
        critical_per_category: Dict of category -> set of critical rsids
        
    Returns:
        Dict of category -> CategoryCoverage
    """
    critical_per_category = critical_per_category or {}
    results = {}
    
    for category, markers in category_markers.items():
        found = 0
        called = 0
        critical_missing = []
        
        critical_set = critical_per_category.get(category, set())
        
        for rsid in markers:
            if rsid in genotypes:
                found += 1
                geno = genotypes[rsid]
                if geno and geno not in ("--", "00", "??", "NC", ""):
                    called += 1
            elif rsid in critical_set:
                critical_missing.append(rsid)
        
        total = len(markers)
        coverage_rate = called / total if total > 0 else 0
        
        # Determine confidence level
        if coverage_rate >= 0.9:
            conf = ConfidenceLevel.DEFINITIVE
        elif coverage_rate >= 0.7:
            conf = ConfidenceLevel.HIGH
        elif coverage_rate >= 0.5:
            conf = ConfidenceLevel.MEDIUM
        elif coverage_rate >= 0.3:
            conf = ConfidenceLevel.LOW
        else:
            conf = ConfidenceLevel.UNCERTAIN
        
        # Downgrade if critical markers missing
        if critical_missing:
            if conf == ConfidenceLevel.DEFINITIVE:
                conf = ConfidenceLevel.HIGH
            elif conf == ConfidenceLevel.HIGH:
                conf = ConfidenceLevel.MEDIUM
        
        results[category] = CategoryCoverage(
            category=category,
            markers_in_reference=total,
            markers_found=found,
            markers_called=called,
            coverage_rate=coverage_rate,
            critical_missing=critical_missing,
            confidence_level=conf
        )
    
    return results


def detect_platform(
    genotypes: Dict[str, str],
    total_snps: int
) -> PlatformDetection:
    """
    Detect genotyping platform from marker profile.
    
    Args:
        genotypes: Dict of rsid -> genotype
        total_snps: Total number of SNPs in file
        
    Returns:
        PlatformDetection object
    """
    best_match = "unknown"
    best_confidence = 0.0
    best_details = None
    alternatives = []
    sig_found = 0
    sig_total = 0
    
    for platform, sig in PLATFORM_SIGNATURES.items():
        min_snps, max_snps = sig["marker_count_range"]
        signature_snps = sig.get("signature_snps", [])
        
        # Check SNP count match
        count_match = 0.0
        if min_snps <= total_snps <= max_snps:
            count_match = 1.0
        elif total_snps < min_snps:
            count_match = total_snps / min_snps
        else:  # total_snps > max_snps
            count_match = max_snps / total_snps
        
        # Check signature SNPs
        sig_matches = sum(1 for snp in signature_snps if snp in genotypes)
        sig_match_rate = sig_matches / len(signature_snps) if signature_snps else 0.5
        
        # Combined confidence
        confidence = 0.6 * count_match + 0.4 * sig_match_rate
        
        if confidence > best_confidence:
            best_confidence = confidence
            best_match = platform
            best_details = sig["description"]
            sig_found = sig_matches
            sig_total = len(signature_snps)
        elif confidence > 0.5:
            alternatives.append(platform)
    
    # Confidence level
    if best_confidence < 0.3:
        best_match = "unknown"
        best_details = "Unable to determine platform"
    
    return PlatformDetection(
        detected_platform=best_match,
        confidence=best_confidence,
        total_snps=total_snps,
        platform_details=best_details,
        alternative_matches=alternatives[:3],
        signature_snps_found=sig_found,
        signature_snps_total=sig_total
    )


def assess_critical_markers(
    genotypes: Dict[str, str]
) -> Tuple[Dict[str, str], List[str]]:
    """
    Check presence and call status of critical markers.
    
    Args:
        genotypes: Dict of rsid -> genotype
        
    Returns:
        Tuple of (status_dict, warning_messages)
    """
    status = {}
    warnings = []
    
    for rsid, info in CRITICAL_MARKERS.items():
        if rsid not in genotypes:
            status[rsid] = "missing"
            warnings.append(f"⚠️ Missing {info['gene']} ({rsid}): {info['reason']}")
        else:
            geno = genotypes[rsid]
            if geno in ("--", "00", "??", "NC", ""):
                status[rsid] = "no_call"
                warnings.append(f"⚠️ No-call {info['gene']} ({rsid}): {info['reason']}")
            else:
                status[rsid] = "present"
    
    return status, warnings


def calculate_overall_score(
    call_rate: float,
    critical_present: int,
    critical_total: int,
    category_coverage_avg: float
) -> Tuple[float, QualityGrade]:
    """
    Calculate overall quality score and grade.
    
    Args:
        call_rate: Overall call rate (0-1)
        critical_present: Number of critical markers present
        critical_total: Total critical markers
        category_coverage_avg: Average category coverage (0-1)
        
    Returns:
        Tuple of (score 0-100, QualityGrade)
    """
    # Weighted scoring
    call_rate_score = call_rate * 40  # 40 points max
    critical_score = (critical_present / critical_total * 30) if critical_total > 0 else 30  # 30 points max
    coverage_score = category_coverage_avg * 30  # 30 points max
    
    total_score = call_rate_score + critical_score + coverage_score
    
    # Determine grade
    if total_score >= 95 and call_rate >= 0.95:
        grade = QualityGrade.EXCELLENT
    elif total_score >= 85 and call_rate >= 0.90:
        grade = QualityGrade.GOOD
    elif total_score >= 70 and call_rate >= 0.80:
        grade = QualityGrade.FAIR
    elif total_score >= 50:
        grade = QualityGrade.POOR
    else:
        grade = QualityGrade.FAIL
    
    return total_score, grade


def generate_quality_report(
    genotypes: Dict[str, str],
    chromosome_map: Dict[str, str] = None,
    category_markers: Dict[str, Set[str]] = None
) -> QualityReport:
    """
    Generate comprehensive quality assessment report.
    
    Args:
        genotypes: Dict of rsid -> genotype
        chromosome_map: Optional dict of rsid -> chromosome
        category_markers: Optional dict of category -> marker rsids
        
    Returns:
        QualityReport object
    """
    total_snps = len(genotypes)
    
    # Count called SNPs
    called_snps = sum(
        1 for geno in genotypes.values()
        if geno and geno not in ("--", "00", "??", "NC", "")
    )
    
    call_rate = called_snps / total_snps if total_snps > 0 else 0
    
    # Chromosome quality (if map provided)
    chrom_quality = {}
    if chromosome_map:
        chrom_quality = assess_chromosome_quality(genotypes, chromosome_map)
    
    # Category coverage (if markers provided)
    cat_coverage = {}
    coverage_avg = 0.8  # Default assumption
    if category_markers:
        cat_coverage = assess_category_coverage(genotypes, category_markers)
        if cat_coverage:
            coverage_avg = sum(c.coverage_rate for c in cat_coverage.values()) / len(cat_coverage)
    
    # Platform detection
    platform = detect_platform(genotypes, total_snps)
    
    # Critical markers
    critical_status, critical_warnings = assess_critical_markers(genotypes)
    critical_present = sum(1 for s in critical_status.values() if s == "present")
    critical_total = len(CRITICAL_MARKERS)
    
    # Overall score
    score, grade = calculate_overall_score(call_rate, critical_present, critical_total, coverage_avg)
    
    # Build warnings and recommendations
    warnings = critical_warnings.copy()
    recommendations = []
    
    if call_rate < 0.95:
        warnings.append(f"Call rate ({call_rate*100:.1f}%) is below optimal (>95%)")
        recommendations.append("Consider re-running or getting new sample if low call rate affects results")
    
    if platform.detected_platform == "unknown":
        warnings.append("Could not identify genotyping platform")
        recommendations.append("Verify file format and source")
    
    missing_critical = [k for k, v in critical_status.items() if v != "present"]
    if missing_critical:
        recommendations.append(f"Consider targeted testing for {len(missing_critical)} missing critical markers")
    
    # Confidence level
    if grade == QualityGrade.EXCELLENT:
        conf = ConfidenceLevel.HIGH
    elif grade == QualityGrade.GOOD:
        conf = ConfidenceLevel.HIGH
    elif grade == QualityGrade.FAIR:
        conf = ConfidenceLevel.MEDIUM
    elif grade == QualityGrade.POOR:
        conf = ConfidenceLevel.LOW
    else:
        conf = ConfidenceLevel.UNCERTAIN
    
    return QualityReport(
        overall_grade=grade,
        overall_score=score,
        total_snps=total_snps,
        called_snps=called_snps,
        overall_call_rate=call_rate,
        chromosome_quality=chrom_quality,
        category_coverage=cat_coverage,
        platform=platform,
        critical_markers_status=critical_status,
        warnings=warnings,
        recommendations=recommendations,
        confidence_level=conf
    )


# =============================================================================
# CONFIDENCE ADJUSTMENTS
# =============================================================================

def adjust_confidence_for_quality(
    base_confidence: ConfidenceLevel,
    quality_report: QualityReport,
    category: str = None
) -> ConfidenceLevel:
    """
    Adjust confidence level based on data quality.
    
    Args:
        base_confidence: Base confidence from analysis
        quality_report: QualityReport from quality assessment
        category: Optional category to check specific coverage
        
    Returns:
        Adjusted ConfidenceLevel
    """
    # Get category-specific confidence if available
    if category and category in quality_report.category_coverage:
        cat_conf = quality_report.category_coverage[category].confidence_level
        # Take minimum of base and category confidence
        levels = [ConfidenceLevel.UNCERTAIN, ConfidenceLevel.LOW, 
                 ConfidenceLevel.MEDIUM, ConfidenceLevel.HIGH, ConfidenceLevel.DEFINITIVE]
        return levels[min(levels.index(base_confidence), levels.index(cat_conf))]
    
    # Adjust based on overall quality
    if quality_report.overall_grade == QualityGrade.FAIL:
        return ConfidenceLevel.UNCERTAIN
    elif quality_report.overall_grade == QualityGrade.POOR:
        if base_confidence in (ConfidenceLevel.HIGH, ConfidenceLevel.DEFINITIVE):
            return ConfidenceLevel.MEDIUM
    elif quality_report.overall_grade == QualityGrade.FAIR:
        if base_confidence == ConfidenceLevel.DEFINITIVE:
            return ConfidenceLevel.HIGH
    
    return base_confidence


def get_quality_summary_text(quality_report: QualityReport) -> str:
    """
    Generate human-readable quality summary.
    
    Args:
        quality_report: QualityReport object
        
    Returns:
        Formatted text summary
    """
    lines = [
        f"📊 DATA QUALITY: Grade {quality_report.overall_grade.value} ({quality_report.overall_score:.0f}/100)",
        f"   Platform: {quality_report.platform.platform_details or quality_report.platform.detected_platform}",
        f"   Total SNPs: {quality_report.total_snps:,}",
        f"   Call Rate: {quality_report.overall_call_rate*100:.1f}%",
        ""
    ]
    
    if quality_report.warnings:
        lines.append("⚠️  WARNINGS:")
        for w in quality_report.warnings[:5]:
            lines.append(f"   • {w}")
        lines.append("")
    
    if quality_report.recommendations:
        lines.append("💡 RECOMMENDATIONS:")
        for r in quality_report.recommendations[:3]:
            lines.append(f"   • {r}")
    
    return "\n".join(lines)
