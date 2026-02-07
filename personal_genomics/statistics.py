"""
Statistical Functions for Personal Genomics Analysis

Provides comprehensive statistical tools for confidence estimation,
interval calculation, and uncertainty quantification.

Author: OpenClaw AI
Date: 2026-02-07
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import (
    Callable, Generic, List, Optional, Sequence, Tuple, TypeVar, Union
)

import numpy as np
from scipy import stats
from scipy.special import betainc, betaln


# =============================================================================
# TYPE DEFINITIONS
# =============================================================================

T = TypeVar('T', bound=float)


class ConfidenceLevel(Enum):
    """Confidence level categories for genomic findings."""
    DEFINITIVE = "DEFINITIVE"  # Defining SNP directly observed
    HIGH = "HIGH"              # Multiple consistent markers
    MEDIUM = "MEDIUM"          # Inferred from related markers
    LOW = "LOW"                # Single marker or chip limitation
    UNCERTAIN = "UNCERTAIN"    # Insufficient data


@dataclass(frozen=True)
class ConfidenceInterval:
    """Immutable confidence interval with metadata."""
    lower: float
    upper: float
    estimate: float
    confidence: float = 0.95
    method: str = "normal"
    n: Optional[int] = None
    
    @property
    def width(self) -> float:
        """Width of the confidence interval."""
        return self.upper - self.lower
    
    @property
    def margin_of_error(self) -> float:
        """Margin of error (half-width)."""
        return self.width / 2
    
    def contains(self, value: float) -> bool:
        """Check if value falls within the interval."""
        return self.lower <= value <= self.upper
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "estimate": round(self.estimate, 4),
            "ci_lower": round(self.lower, 4),
            "ci_upper": round(self.upper, 4),
            "confidence": self.confidence,
            "method": self.method,
            "n": self.n
        }


@dataclass
class StatisticalResult:
    """Result container with full statistical metadata."""
    value: float
    ci_lower: float
    ci_upper: float
    confidence_level: ConfidenceLevel
    p_value: Optional[float] = None
    standard_error: Optional[float] = None
    n_observations: Optional[int] = None
    n_total: Optional[int] = None
    method: str = "standard"
    quality_score: Optional[float] = None
    notes: List[str] = field(default_factory=list)
    
    @property
    def confidence_numeric(self) -> float:
        """Convert confidence level to numeric score 0-1."""
        mapping = {
            ConfidenceLevel.DEFINITIVE: 0.99,
            ConfidenceLevel.HIGH: 0.90,
            ConfidenceLevel.MEDIUM: 0.70,
            ConfidenceLevel.LOW: 0.50,
            ConfidenceLevel.UNCERTAIN: 0.30,
        }
        return mapping.get(self.confidence_level, 0.50)
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        result = {
            "value": round(self.value, 4),
            "ci_lower": round(self.ci_lower, 4),
            "ci_upper": round(self.ci_upper, 4),
            "confidence": self.confidence_level.value,
            "confidence_numeric": self.confidence_numeric,
        }
        if self.p_value is not None:
            result["p_value"] = round(self.p_value, 6)
        if self.standard_error is not None:
            result["standard_error"] = round(self.standard_error, 4)
        if self.n_observations is not None:
            result["n_markers"] = self.n_observations
        if self.n_total is not None:
            result["n_total"] = self.n_total
        if self.quality_score is not None:
            result["quality_score"] = round(self.quality_score, 3)
        if self.notes:
            result["notes"] = self.notes
        return result


# =============================================================================
# CORE STATISTICAL FUNCTIONS
# =============================================================================

def confidence_interval(
    estimate: float,
    n: int,
    confidence: float = 0.95,
    std: Optional[float] = None,
    method: str = "normal"
) -> ConfidenceInterval:
    """
    Calculate confidence interval for a point estimate.
    
    Args:
        estimate: Point estimate (e.g., sample mean, proportion)
        n: Sample size
        confidence: Confidence level (default 0.95)
        std: Standard deviation (estimated if not provided)
        method: "normal" for large samples, "t" for small samples
        
    Returns:
        ConfidenceInterval object
        
    Example:
        >>> ci = confidence_interval(0.65, n=100, confidence=0.95)
        >>> print(f"{ci.lower:.3f} - {ci.upper:.3f}")
        0.557 - 0.743
    """
    if n <= 0:
        return ConfidenceInterval(
            lower=estimate, upper=estimate, estimate=estimate,
            confidence=confidence, method="insufficient_data", n=0
        )
    
    # Estimate standard error
    if std is not None:
        se = std / math.sqrt(n)
    else:
        # For proportions, use sqrt(p*(1-p)/n)
        if 0 <= estimate <= 1:
            se = math.sqrt(estimate * (1 - estimate) / n) if n > 0 else 0
        else:
            # Fallback: assume some variance
            se = abs(estimate) * 0.1 / math.sqrt(n) if n > 0 else 0
    
    # Calculate critical value
    alpha = 1 - confidence
    if method == "t" and n < 30:
        critical = stats.t.ppf(1 - alpha / 2, df=n - 1)
    else:
        critical = stats.norm.ppf(1 - alpha / 2)
    
    margin = critical * se
    
    return ConfidenceInterval(
        lower=estimate - margin,
        upper=estimate + margin,
        estimate=estimate,
        confidence=confidence,
        method=method,
        n=n
    )


def wilson_score_interval(
    successes: int,
    n: int,
    confidence: float = 0.95
) -> ConfidenceInterval:
    """
    Wilson score interval for proportions.
    
    More accurate than normal approximation, especially for small samples
    or proportions near 0 or 1.
    
    Args:
        successes: Number of successes
        n: Total trials
        confidence: Confidence level (default 0.95)
        
    Returns:
        ConfidenceInterval object
        
    Example:
        >>> ci = wilson_score_interval(8, 10)  # 80% success
        >>> print(f"{ci.lower:.3f} - {ci.upper:.3f}")
        0.493 - 0.943
    """
    if n <= 0:
        return ConfidenceInterval(
            lower=0.0, upper=1.0, estimate=0.5,
            confidence=confidence, method="wilson_empty", n=0
        )
    
    p_hat = successes / n
    z = stats.norm.ppf(1 - (1 - confidence) / 2)
    z2 = z ** 2
    
    denominator = 1 + z2 / n
    center = (p_hat + z2 / (2 * n)) / denominator
    margin = (z / denominator) * math.sqrt(
        p_hat * (1 - p_hat) / n + z2 / (4 * n ** 2)
    )
    
    lower = max(0.0, center - margin)
    upper = min(1.0, center + margin)
    
    return ConfidenceInterval(
        lower=lower,
        upper=upper,
        estimate=p_hat,
        confidence=confidence,
        method="wilson",
        n=n
    )


def bayesian_posterior(
    prior_alpha: float,
    prior_beta: float,
    successes: int,
    failures: int
) -> Tuple[float, ConfidenceInterval]:
    """
    Bayesian posterior for beta-binomial model.
    
    Args:
        prior_alpha: Prior alpha parameter (e.g., 1 for uniform)
        prior_beta: Prior beta parameter (e.g., 1 for uniform)
        successes: Observed successes
        failures: Observed failures
        
    Returns:
        Tuple of (posterior mean, 95% credible interval)
        
    Example:
        >>> mean, ci = bayesian_posterior(1, 1, 8, 2)  # Uniform prior, 8/10 success
        >>> print(f"Posterior mean: {mean:.3f}, CI: [{ci.lower:.3f}, {ci.upper:.3f}]")
        Posterior mean: 0.750, CI: [0.473, 0.927]
    """
    post_alpha = prior_alpha + successes
    post_beta = prior_beta + failures
    
    # Posterior mean
    mean = post_alpha / (post_alpha + post_beta)
    
    # 95% credible interval using beta quantiles
    lower = stats.beta.ppf(0.025, post_alpha, post_beta)
    upper = stats.beta.ppf(0.975, post_alpha, post_beta)
    
    ci = ConfidenceInterval(
        lower=lower,
        upper=upper,
        estimate=mean,
        confidence=0.95,
        method="bayesian_beta",
        n=successes + failures
    )
    
    return mean, ci


def bootstrap_ci(
    data: Sequence[float],
    statistic: Callable[[Sequence[float]], float],
    n_bootstrap: int = 1000,
    confidence: float = 0.95,
    random_state: Optional[int] = None
) -> ConfidenceInterval:
    """
    Bootstrap confidence interval for any statistic.
    
    Args:
        data: Sequence of observations
        statistic: Function to compute statistic from sample
        n_bootstrap: Number of bootstrap iterations (default 1000)
        confidence: Confidence level (default 0.95)
        random_state: Random seed for reproducibility
        
    Returns:
        ConfidenceInterval object
        
    Example:
        >>> data = [0.1, 0.2, 0.15, 0.3, 0.25, 0.18]
        >>> ci = bootstrap_ci(data, np.median, n_bootstrap=1000)
        >>> print(f"Median: {ci.estimate:.3f}, CI: [{ci.lower:.3f}, {ci.upper:.3f}]")
    """
    if len(data) == 0:
        return ConfidenceInterval(
            lower=0.0, upper=0.0, estimate=0.0,
            confidence=confidence, method="bootstrap_empty", n=0
        )
    
    rng = np.random.default_rng(random_state)
    data_arr = np.array(data)
    n = len(data_arr)
    
    # Original statistic
    original = statistic(data_arr)
    
    # Bootstrap resampling
    bootstrap_stats = []
    for _ in range(n_bootstrap):
        sample = rng.choice(data_arr, size=n, replace=True)
        bootstrap_stats.append(statistic(sample))
    
    bootstrap_stats = np.array(bootstrap_stats)
    
    # Percentile method
    alpha = 1 - confidence
    lower = np.percentile(bootstrap_stats, 100 * alpha / 2)
    upper = np.percentile(bootstrap_stats, 100 * (1 - alpha / 2))
    
    return ConfidenceInterval(
        lower=float(lower),
        upper=float(upper),
        estimate=float(original),
        confidence=confidence,
        method="bootstrap_percentile",
        n=n
    )


def effect_size_ci(
    beta: float,
    se: float,
    confidence: float = 0.95
) -> ConfidenceInterval:
    """
    Confidence interval for effect size (e.g., GWAS beta).
    
    Args:
        beta: Effect size estimate (log odds ratio or beta coefficient)
        se: Standard error of the estimate
        confidence: Confidence level (default 0.95)
        
    Returns:
        ConfidenceInterval object
        
    Example:
        >>> ci = effect_size_ci(0.25, 0.08)  # Beta=0.25, SE=0.08
        >>> print(f"Effect: {ci.estimate:.3f}, CI: [{ci.lower:.3f}, {ci.upper:.3f}]")
        Effect: 0.250, CI: [0.093, 0.407]
    """
    z = stats.norm.ppf(1 - (1 - confidence) / 2)
    margin = z * se
    
    return ConfidenceInterval(
        lower=beta - margin,
        upper=beta + margin,
        estimate=beta,
        confidence=confidence,
        method="effect_size",
        n=None
    )


def marker_coverage_weight(
    found: int,
    total: int,
    min_threshold: float = 0.3,
    optimal_threshold: float = 0.7
) -> Tuple[float, ConfidenceLevel]:
    """
    Calculate quality score and confidence level based on marker coverage.
    
    Args:
        found: Number of markers found in data
        total: Total markers in reference panel
        min_threshold: Minimum coverage for valid result (default 0.3)
        optimal_threshold: Coverage for high confidence (default 0.7)
        
    Returns:
        Tuple of (quality_score, confidence_level)
        
    Example:
        >>> score, level = marker_coverage_weight(15, 20)
        >>> print(f"Quality: {score:.2f}, Confidence: {level.value}")
        Quality: 0.75, Confidence: HIGH
    """
    if total == 0:
        return 0.0, ConfidenceLevel.UNCERTAIN
    
    coverage = found / total
    
    # Quality score is the coverage ratio
    quality_score = coverage
    
    # Determine confidence level
    if coverage >= 0.9:
        level = ConfidenceLevel.DEFINITIVE
    elif coverage >= optimal_threshold:
        level = ConfidenceLevel.HIGH
    elif coverage >= 0.5:
        level = ConfidenceLevel.MEDIUM
    elif coverage >= min_threshold:
        level = ConfidenceLevel.LOW
    else:
        level = ConfidenceLevel.UNCERTAIN
    
    return quality_score, level


# =============================================================================
# P-VALUE CALCULATIONS
# =============================================================================

def proportion_test_pvalue(
    observed: float,
    expected: float,
    n: int,
    alternative: str = "two-sided"
) -> float:
    """
    Test if observed proportion differs from expected.
    
    Args:
        observed: Observed proportion
        expected: Expected proportion under null hypothesis
        n: Sample size
        alternative: "two-sided", "greater", or "less"
        
    Returns:
        p-value
        
    Example:
        >>> p = proportion_test_pvalue(0.65, 0.5, n=100)
        >>> print(f"p-value: {p:.4f}")  # Is 65% significantly different from 50%?
        p-value: 0.0027
    """
    if n <= 0 or expected <= 0 or expected >= 1:
        return 1.0
    
    # Standard error under null
    se = math.sqrt(expected * (1 - expected) / n)
    
    if se == 0:
        return 1.0 if observed == expected else 0.0
    
    z = (observed - expected) / se
    
    if alternative == "two-sided":
        return 2 * (1 - stats.norm.cdf(abs(z)))
    elif alternative == "greater":
        return 1 - stats.norm.cdf(z)
    else:  # less
        return stats.norm.cdf(z)


def chi2_test_pvalue(observed: Sequence[int], expected: Sequence[float]) -> float:
    """
    Chi-squared goodness of fit test.
    
    Args:
        observed: Observed counts
        expected: Expected counts
        
    Returns:
        p-value
        
    Example:
        >>> observed = [45, 35, 20]  # Observed genotype counts
        >>> expected = [50, 30, 20]  # Expected under HWE
        >>> p = chi2_test_pvalue(observed, expected)
    """
    obs = np.array(observed)
    exp = np.array(expected)
    
    if len(obs) != len(exp) or len(obs) < 2:
        return 1.0
    
    # Avoid division by zero
    exp = np.maximum(exp, 1e-10)
    
    chi2_stat = np.sum((obs - exp) ** 2 / exp)
    df = len(obs) - 1
    
    return 1 - stats.chi2.cdf(chi2_stat, df)


def binomial_test_pvalue(
    successes: int,
    n: int,
    p: float = 0.5,
    alternative: str = "two-sided"
) -> float:
    """
    Exact binomial test.
    
    Args:
        successes: Number of successes
        n: Total trials
        p: Expected probability under null
        alternative: "two-sided", "greater", or "less"
        
    Returns:
        p-value
        
    Example:
        >>> p = binomial_test_pvalue(8, 10, p=0.5)  # 8/10 vs 50%
        >>> print(f"p-value: {p:.4f}")
    """
    result = stats.binomtest(successes, n, p, alternative=alternative)
    return result.pvalue


# =============================================================================
# GENOMICS-SPECIFIC STATISTICS
# =============================================================================

def prs_percentile_ci(
    raw_score: float,
    n_markers: int,
    total_markers: int,
    population_mean: float = 0.0,
    population_sd: float = 1.0,
    confidence: float = 0.95
) -> StatisticalResult:
    """
    Calculate PRS percentile with confidence interval.
    
    Args:
        raw_score: Raw polygenic risk score
        n_markers: Number of markers used
        total_markers: Total markers in full PRS
        population_mean: Population mean (default 0)
        population_sd: Population SD (default 1)
        confidence: Confidence level
        
    Returns:
        StatisticalResult with percentile and CI
        
    Example:
        >>> result = prs_percentile_ci(1.5, n_markers=50, total_markers=100)
        >>> print(f"Percentile: {result.value:.1f}, CI: [{result.ci_lower:.1f}, {result.ci_upper:.1f}]")
    """
    if total_markers == 0 or n_markers == 0:
        return StatisticalResult(
            value=50.0,
            ci_lower=0.0,
            ci_upper=100.0,
            confidence_level=ConfidenceLevel.UNCERTAIN,
            n_observations=0,
            n_total=total_markers,
            notes=["Insufficient markers for calculation"]
        )
    
    coverage = n_markers / total_markers
    
    # Standardize
    z_score = (raw_score - population_mean) / population_sd if population_sd > 0 else 0
    
    # Convert to percentile
    percentile = stats.norm.cdf(z_score) * 100
    
    # Estimate standard error - increases with lower coverage
    # SE scales roughly as sqrt(total/found) due to sampling variance
    se_multiplier = math.sqrt(total_markers / n_markers) if n_markers > 0 else 3.0
    se_percentile = 10 * se_multiplier  # Base SE of ~10 percentile points
    
    # Calculate CI
    z = stats.norm.ppf(1 - (1 - confidence) / 2)
    lower = max(0, percentile - z * se_percentile)
    upper = min(100, percentile + z * se_percentile)
    
    # Determine confidence level
    quality_score, conf_level = marker_coverage_weight(n_markers, total_markers)
    
    notes = []
    if coverage < 0.5:
        notes.append(f"LOW COVERAGE: Only {coverage*100:.0f}% of markers available")
    
    return StatisticalResult(
        value=percentile,
        ci_lower=lower,
        ci_upper=upper,
        confidence_level=conf_level,
        standard_error=se_percentile,
        n_observations=n_markers,
        n_total=total_markers,
        quality_score=quality_score,
        method="prs_percentile",
        notes=notes
    )


def ancestry_similarity_stats(
    matching_alleles: int,
    total_alleles: int,
    baseline_frequency: float = 0.5
) -> StatisticalResult:
    """
    Calculate ancestry similarity with statistical significance.
    
    Args:
        matching_alleles: Number of alleles matching reference population
        total_alleles: Total alleles compared
        baseline_frequency: Expected frequency under random (typically 0.5)
        
    Returns:
        StatisticalResult with similarity score, CI, and p-value
        
    Example:
        >>> result = ancestry_similarity_stats(16, 24, baseline_frequency=0.5)
        >>> print(f"Similarity: {result.value:.1f}%, p-value: {result.p_value:.4f}")
    """
    if total_alleles == 0:
        return StatisticalResult(
            value=0.0,
            ci_lower=0.0,
            ci_upper=100.0,
            confidence_level=ConfidenceLevel.UNCERTAIN,
            n_observations=0,
            notes=["No markers compared"]
        )
    
    # Calculate proportion
    similarity = matching_alleles / total_alleles
    
    # Wilson score interval (better for proportions)
    ci = wilson_score_interval(matching_alleles, total_alleles)
    
    # P-value: is this different from random?
    p_value = proportion_test_pvalue(similarity, baseline_frequency, total_alleles)
    
    # Standard error
    se = math.sqrt(similarity * (1 - similarity) / total_alleles)
    
    # Confidence level based on markers and significance
    if total_alleles >= 50 and p_value < 0.001:
        conf_level = ConfidenceLevel.HIGH
    elif total_alleles >= 20 and p_value < 0.05:
        conf_level = ConfidenceLevel.MEDIUM
    elif total_alleles >= 10:
        conf_level = ConfidenceLevel.LOW
    else:
        conf_level = ConfidenceLevel.UNCERTAIN
    
    return StatisticalResult(
        value=similarity * 100,  # As percentage
        ci_lower=ci.lower * 100,
        ci_upper=ci.upper * 100,
        confidence_level=conf_level,
        p_value=p_value,
        standard_error=se * 100,
        n_observations=total_alleles,
        method="wilson_score"
    )


def haplogroup_confidence(
    defining_snps_found: int,
    defining_snps_total: int,
    consistent_markers: int,
    inconsistent_markers: int
) -> Tuple[ConfidenceLevel, float]:
    """
    Determine haplogroup call confidence.
    
    Args:
        defining_snps_found: Number of clade-defining SNPs observed
        defining_snps_total: Total defining SNPs for this haplogroup
        consistent_markers: Markers consistent with the call
        inconsistent_markers: Markers inconsistent with the call
        
    Returns:
        Tuple of (confidence_level, confidence_score)
        
    Example:
        >>> level, score = haplogroup_confidence(
        ...     defining_snps_found=1, defining_snps_total=1,
        ...     consistent_markers=5, inconsistent_markers=0
        ... )
        >>> print(f"Confidence: {level.value}, Score: {score:.2f}")
    """
    total_markers = consistent_markers + inconsistent_markers
    
    # If we found the defining SNP
    if defining_snps_found >= 1 and defining_snps_total >= 1:
        defining_ratio = defining_snps_found / defining_snps_total
        
        if defining_ratio >= 1.0 and inconsistent_markers == 0:
            return ConfidenceLevel.DEFINITIVE, 0.99
        elif defining_ratio >= 0.5:
            return ConfidenceLevel.HIGH, 0.85 + 0.10 * defining_ratio
    
    # No defining SNP found - infer from related markers
    if total_markers == 0:
        return ConfidenceLevel.UNCERTAIN, 0.30
    
    consistency_ratio = consistent_markers / total_markers
    
    if consistency_ratio >= 0.9 and consistent_markers >= 3:
        return ConfidenceLevel.HIGH, 0.80 + 0.10 * consistency_ratio
    elif consistency_ratio >= 0.7 and consistent_markers >= 2:
        return ConfidenceLevel.MEDIUM, 0.60 + 0.15 * consistency_ratio
    elif consistent_markers >= 1:
        return ConfidenceLevel.LOW, 0.40 + 0.20 * consistency_ratio
    else:
        return ConfidenceLevel.UNCERTAIN, 0.30


def diplotype_confidence(
    star_allele_1: str,
    star_allele_2: str,
    markers_for_allele_1: int,
    total_markers_allele_1: int,
    markers_for_allele_2: int,
    total_markers_allele_2: int,
    phase_ambiguous: bool = False
) -> Tuple[float, List[str]]:
    """
    Calculate confidence for a diplotype (star allele) call.
    
    Args:
        star_allele_1: First star allele (e.g., "*1")
        star_allele_2: Second star allele (e.g., "*4")
        markers_for_allele_1: Markers found for first allele
        total_markers_allele_1: Total markers defining first allele
        markers_for_allele_2: Markers found for second allele
        total_markers_allele_2: Total markers defining second allele
        phase_ambiguous: Whether phase cannot be determined
        
    Returns:
        Tuple of (confidence_score, warning_messages)
        
    Example:
        >>> conf, warnings = diplotype_confidence("*1", "*4", 3, 3, 2, 3, phase_ambiguous=False)
        >>> print(f"Confidence: {conf:.2f}")
    """
    warnings = []
    
    # Calculate coverage for each allele
    cov_1 = markers_for_allele_1 / total_markers_allele_1 if total_markers_allele_1 > 0 else 0
    cov_2 = markers_for_allele_2 / total_markers_allele_2 if total_markers_allele_2 > 0 else 0
    
    # Base confidence from marker coverage
    base_conf = (cov_1 + cov_2) / 2
    
    # Penalties
    if phase_ambiguous:
        base_conf *= 0.85
        warnings.append("Phase ambiguous - diplotype inferred")
    
    if cov_1 < 0.5:
        base_conf *= 0.9
        warnings.append(f"Low coverage for {star_allele_1} ({cov_1*100:.0f}%)")
    
    if cov_2 < 0.5:
        base_conf *= 0.9
        warnings.append(f"Low coverage for {star_allele_2} ({cov_2*100:.0f}%)")
    
    # Check for reference allele (*1) which is often defined by absence
    if star_allele_1 == "*1" or star_allele_2 == "*1":
        warnings.append("*1 allele inferred by absence of other variants")
        base_conf = min(base_conf, 0.9)  # Cap at 90% for *1 inference
    
    return min(1.0, base_conf), warnings


def trait_probability_ci(
    genotype: str,
    effect_allele: str,
    baseline_probability: float,
    odds_ratio: float,
    odds_ratio_ci_lower: float = None,
    odds_ratio_ci_upper: float = None
) -> StatisticalResult:
    """
    Calculate trait probability with confidence interval from OR.
    
    Args:
        genotype: User's genotype (e.g., "AG")
        effect_allele: Effect allele (e.g., "A")
        baseline_probability: Population baseline probability
        odds_ratio: Odds ratio for effect allele
        odds_ratio_ci_lower: Lower bound of OR CI
        odds_ratio_ci_upper: Upper bound of OR CI
        
    Returns:
        StatisticalResult with probability and CI
        
    Example:
        >>> result = trait_probability_ci("AG", "A", 0.3, 2.5, 1.8, 3.4)
        >>> print(f"Probability: {result.value:.2f}, CI: [{result.ci_lower:.2f}, {result.ci_upper:.2f}]")
    """
    # Count effect alleles
    effect_count = genotype.upper().count(effect_allele.upper())
    
    # Calculate odds for this genotype
    baseline_odds = baseline_probability / (1 - baseline_probability) if baseline_probability < 1 else 100
    
    # Apply OR based on allele count (additive model)
    adjusted_odds = baseline_odds * (odds_ratio ** effect_count)
    
    # Convert back to probability
    probability = adjusted_odds / (1 + adjusted_odds)
    
    # Calculate CI bounds if OR CI provided
    if odds_ratio_ci_lower and odds_ratio_ci_upper:
        lower_odds = baseline_odds * (odds_ratio_ci_lower ** effect_count)
        upper_odds = baseline_odds * (odds_ratio_ci_upper ** effect_count)
        
        ci_lower = lower_odds / (1 + lower_odds)
        ci_upper = upper_odds / (1 + upper_odds)
    else:
        # Estimate CI from effect size uncertainty
        ci = confidence_interval(probability, n=1000)  # Rough estimate
        ci_lower = ci.lower
        ci_upper = ci.upper
    
    # Confidence level
    if odds_ratio_ci_lower and odds_ratio_ci_upper:
        # Check if CI excludes 1.0
        if odds_ratio_ci_lower > 1.0 or odds_ratio_ci_upper < 1.0:
            conf_level = ConfidenceLevel.HIGH
        else:
            conf_level = ConfidenceLevel.MEDIUM
    else:
        conf_level = ConfidenceLevel.LOW
    
    return StatisticalResult(
        value=probability,
        ci_lower=max(0, ci_lower),
        ci_upper=min(1, ci_upper),
        confidence_level=conf_level,
        method="odds_ratio_transform"
    )


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def combine_confidence_levels(levels: Sequence[ConfidenceLevel]) -> ConfidenceLevel:
    """
    Combine multiple confidence levels (conservative - takes minimum).
    
    Args:
        levels: Sequence of ConfidenceLevel values
        
    Returns:
        Combined ConfidenceLevel (minimum of inputs)
    """
    if not levels:
        return ConfidenceLevel.UNCERTAIN
    
    order = [
        ConfidenceLevel.UNCERTAIN,
        ConfidenceLevel.LOW,
        ConfidenceLevel.MEDIUM,
        ConfidenceLevel.HIGH,
        ConfidenceLevel.DEFINITIVE
    ]
    
    min_idx = min(order.index(level) for level in levels)
    return order[min_idx]


def confidence_to_color(level: ConfidenceLevel) -> str:
    """
    Convert confidence level to color code for visualization.
    
    Args:
        level: ConfidenceLevel enum value
        
    Returns:
        Hex color code
    """
    colors = {
        ConfidenceLevel.DEFINITIVE: "#16a34a",  # Green
        ConfidenceLevel.HIGH: "#22c55e",        # Light green
        ConfidenceLevel.MEDIUM: "#eab308",      # Yellow
        ConfidenceLevel.LOW: "#f97316",         # Orange
        ConfidenceLevel.UNCERTAIN: "#dc2626",   # Red
    }
    return colors.get(level, "#6b7280")  # Gray default


def format_ci_string(
    estimate: float,
    lower: float,
    upper: float,
    decimals: int = 1,
    as_percentage: bool = False
) -> str:
    """
    Format confidence interval as string.
    
    Args:
        estimate: Point estimate
        lower: Lower bound
        upper: Upper bound
        decimals: Decimal places
        as_percentage: Whether to format as percentage
        
    Returns:
        Formatted string like "65.0% (58.2-71.8%)"
    """
    if as_percentage:
        return f"{estimate:.{decimals}f}% ({lower:.{decimals}f}-{upper:.{decimals}f}%)"
    else:
        return f"{estimate:.{decimals}f} ({lower:.{decimals}f}-{upper:.{decimals}f})"
