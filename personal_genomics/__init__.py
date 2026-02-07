"""
Personal Genomics Package

Statistical functions and data quality tools for genetic analysis.
"""

from .statistics import (
    # Enums and Types
    ConfidenceLevel,
    ConfidenceInterval,
    StatisticalResult,
    
    # Core statistical functions
    confidence_interval,
    wilson_score_interval,
    bayesian_posterior,
    bootstrap_ci,
    effect_size_ci,
    marker_coverage_weight,
    
    # P-value calculations
    proportion_test_pvalue,
    chi2_test_pvalue,
    binomial_test_pvalue,
    
    # Genomics-specific
    prs_percentile_ci,
    ancestry_similarity_stats,
    haplogroup_confidence,
    diplotype_confidence,
    trait_probability_ci,
    
    # Utilities
    combine_confidence_levels,
    confidence_to_color,
    format_ci_string,
)

from .quality import (
    # Types
    QualityGrade,
    ChromosomeQuality,
    CategoryCoverage,
    PlatformDetection,
    QualityReport,
    
    # Assessment functions
    assess_chromosome_quality,
    assess_category_coverage,
    detect_platform,
    assess_critical_markers,
    calculate_overall_score,
    generate_quality_report,
    adjust_confidence_for_quality,
    get_quality_summary_text,
    
    # Constants
    CRITICAL_MARKERS,
    CHROMOSOME_EXPECTED,
    PLATFORM_SIGNATURES,
)

__all__ = [
    # Statistics
    "ConfidenceLevel",
    "ConfidenceInterval", 
    "StatisticalResult",
    "confidence_interval",
    "wilson_score_interval",
    "bayesian_posterior",
    "bootstrap_ci",
    "effect_size_ci",
    "marker_coverage_weight",
    "proportion_test_pvalue",
    "chi2_test_pvalue",
    "binomial_test_pvalue",
    "prs_percentile_ci",
    "ancestry_similarity_stats",
    "haplogroup_confidence",
    "diplotype_confidence",
    "trait_probability_ci",
    "combine_confidence_levels",
    "confidence_to_color",
    "format_ci_string",
    
    # Quality
    "QualityGrade",
    "ChromosomeQuality",
    "CategoryCoverage",
    "PlatformDetection",
    "QualityReport",
    "assess_chromosome_quality",
    "assess_category_coverage",
    "detect_platform",
    "assess_critical_markers",
    "calculate_overall_score",
    "generate_quality_report",
    "adjust_confidence_for_quality",
    "get_quality_summary_text",
    "CRITICAL_MARKERS",
    "CHROMOSOME_EXPECTED",
    "PLATFORM_SIGNATURES",
]
