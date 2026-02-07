"""
Comprehensive Tests for Statistical Functions

Tests all statistical functions in personal_genomics.statistics module
including confidence intervals, p-values, and genomics-specific calculations.

Author: OpenClaw AI
Date: 2026-02-07
"""

import math
import pytest
import sys
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from personal_genomics.statistics import (
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


# =============================================================================
# CONFIDENCE INTERVAL TESTS
# =============================================================================

class TestConfidenceInterval:
    """Tests for confidence_interval function."""
    
    def test_basic_ci(self):
        """Test basic confidence interval calculation."""
        ci = confidence_interval(0.65, n=100, confidence=0.95)
        
        assert ci.estimate == 0.65
        assert 0.55 < ci.lower < 0.65
        assert 0.65 < ci.upper < 0.75
        assert ci.confidence == 0.95
        assert ci.n == 100
    
    def test_ci_width_decreases_with_n(self):
        """Larger sample size should give narrower CI."""
        ci_small = confidence_interval(0.5, n=10)
        ci_large = confidence_interval(0.5, n=1000)
        
        assert ci_large.width < ci_small.width
    
    def test_ci_contains_estimate(self):
        """CI should always contain the point estimate."""
        ci = confidence_interval(0.75, n=50)
        
        assert ci.contains(0.75)
    
    def test_ci_to_dict(self):
        """Test dictionary conversion."""
        ci = confidence_interval(0.5, n=100)
        d = ci.to_dict()
        
        assert "estimate" in d
        assert "ci_lower" in d
        assert "ci_upper" in d
        assert d["confidence"] == 0.95
    
    def test_empty_sample(self):
        """Test with n=0."""
        ci = confidence_interval(0.5, n=0)
        
        assert ci.n == 0
        assert ci.method == "insufficient_data"


class TestWilsonScoreInterval:
    """Tests for Wilson score interval (for proportions)."""
    
    def test_basic_wilson(self):
        """Test basic Wilson score calculation."""
        ci = wilson_score_interval(8, 10)  # 80% success
        
        assert ci.estimate == 0.8
        assert 0.4 < ci.lower < 0.8  # Wilson is asymmetric
        assert 0.8 < ci.upper < 1.0
        assert ci.method == "wilson"
    
    def test_wilson_zero_successes(self):
        """Test with 0 successes."""
        ci = wilson_score_interval(0, 10)
        
        assert ci.estimate == 0.0
        assert ci.lower == 0.0
        assert ci.upper > 0  # Wilson gives non-zero upper bound
    
    def test_wilson_all_successes(self):
        """Test with 100% success rate."""
        ci = wilson_score_interval(10, 10)
        
        assert ci.estimate == 1.0
        assert ci.lower < 1.0  # Wilson gives less than 100%
        assert ci.upper >= 0.999  # Wilson upper bound capped near 1.0
    
    def test_wilson_empty(self):
        """Test with n=0."""
        ci = wilson_score_interval(0, 0)
        
        assert ci.n == 0


class TestBayesianPosterior:
    """Tests for Bayesian posterior calculation."""
    
    def test_uniform_prior(self):
        """Test with uniform prior (Beta(1,1))."""
        mean, ci = bayesian_posterior(1, 1, 8, 2)  # 8 successes, 2 failures
        
        assert 0.7 < mean < 0.8  # Posterior mean
        assert ci.lower < mean < ci.upper
    
    def test_informative_prior(self):
        """Test with informative prior."""
        mean1, ci1 = bayesian_posterior(1, 1, 5, 5)  # Uninformative
        mean2, ci2 = bayesian_posterior(10, 10, 5, 5)  # Prior pulls toward 0.5
        
        # Informative prior should have narrower CI
        assert ci2.width < ci1.width
    
    def test_posterior_mean_formula(self):
        """Verify posterior mean follows Beta conjugate formula."""
        mean, _ = bayesian_posterior(2, 3, 4, 5)
        expected = (2 + 4) / (2 + 3 + 4 + 5)
        
        assert abs(mean - expected) < 0.001


class TestBootstrapCI:
    """Tests for bootstrap confidence interval."""
    
    def test_bootstrap_mean(self):
        """Test bootstrap CI for mean."""
        import numpy as np
        
        data = [10, 12, 15, 11, 14, 13, 16, 12, 14, 15]
        ci = bootstrap_ci(data, np.mean, n_bootstrap=500, random_state=42)
        
        assert 11 < ci.estimate < 15
        assert ci.lower < ci.estimate < ci.upper
    
    def test_bootstrap_median(self):
        """Test bootstrap CI for median."""
        import numpy as np
        
        data = [1, 2, 3, 100]  # Skewed data
        ci = bootstrap_ci(data, np.median, n_bootstrap=500)
        
        assert ci.estimate == np.median(data)
    
    def test_bootstrap_empty(self):
        """Test with empty data."""
        ci = bootstrap_ci([], lambda x: 0)
        
        assert ci.n == 0


class TestEffectSizeCI:
    """Tests for effect size confidence interval."""
    
    def test_basic_effect_ci(self):
        """Test basic effect size CI."""
        ci = effect_size_ci(0.25, se=0.08)
        
        assert ci.estimate == 0.25
        assert 0.09 < ci.lower < 0.25
        assert 0.25 < ci.upper < 0.41
    
    def test_effect_ci_significance(self):
        """CI not containing 0 indicates significance."""
        ci_sig = effect_size_ci(0.30, se=0.10)
        ci_nonsig = effect_size_ci(0.10, se=0.10)
        
        assert not ci_sig.contains(0)  # Significant
        assert ci_nonsig.contains(0)   # Not significant


# =============================================================================
# MARKER COVERAGE TESTS
# =============================================================================

class TestMarkerCoverageWeight:
    """Tests for marker coverage quality scoring."""
    
    def test_high_coverage(self):
        """Test high marker coverage."""
        score, level = marker_coverage_weight(18, 20)
        
        assert score == 0.9
        assert level == ConfidenceLevel.DEFINITIVE
    
    def test_good_coverage(self):
        """Test good marker coverage."""
        score, level = marker_coverage_weight(15, 20)
        
        assert score == 0.75
        assert level == ConfidenceLevel.HIGH
    
    def test_medium_coverage(self):
        """Test medium coverage."""
        score, level = marker_coverage_weight(10, 20)
        
        assert score == 0.5
        assert level == ConfidenceLevel.MEDIUM
    
    def test_low_coverage(self):
        """Test low coverage."""
        score, level = marker_coverage_weight(6, 20)
        
        assert score == 0.3
        assert level == ConfidenceLevel.LOW
    
    def test_zero_total(self):
        """Test with zero total markers."""
        score, level = marker_coverage_weight(0, 0)
        
        assert score == 0.0
        assert level == ConfidenceLevel.UNCERTAIN


# =============================================================================
# P-VALUE TESTS
# =============================================================================

class TestProportionTestPvalue:
    """Tests for proportion test p-value."""
    
    def test_significant_difference(self):
        """Test clearly significant difference."""
        p = proportion_test_pvalue(0.70, expected=0.50, n=100)
        
        assert p < 0.01  # Should be highly significant
    
    def test_no_difference(self):
        """Test when observed equals expected."""
        p = proportion_test_pvalue(0.50, expected=0.50, n=100)
        
        assert p == 1.0
    
    def test_one_sided_greater(self):
        """Test one-sided test (greater)."""
        p = proportion_test_pvalue(0.60, expected=0.50, n=100, alternative="greater")
        
        assert p < proportion_test_pvalue(0.60, expected=0.50, n=100, alternative="two-sided")
    
    def test_small_sample(self):
        """Test with small sample size."""
        p = proportion_test_pvalue(0.70, expected=0.50, n=10)
        
        # Larger p-value with smaller n
        assert p > proportion_test_pvalue(0.70, expected=0.50, n=100)


class TestBinomialTestPvalue:
    """Tests for exact binomial test."""
    
    def test_coin_flip(self):
        """Test fair coin hypothesis."""
        p = binomial_test_pvalue(8, 10, p=0.5)
        
        assert 0.05 < p < 0.20  # Not highly significant with n=10
    
    def test_extreme_result(self):
        """Test extreme result."""
        p = binomial_test_pvalue(10, 10, p=0.5)
        
        assert p < 0.01  # Very unlikely under fair coin


# =============================================================================
# GENOMICS-SPECIFIC TESTS
# =============================================================================

class TestPRSPercentileCI:
    """Tests for PRS percentile confidence interval."""
    
    def test_full_coverage(self):
        """Test with full marker coverage."""
        result = prs_percentile_ci(1.5, n_markers=100, total_markers=100)
        
        assert result.confidence_level in (ConfidenceLevel.HIGH, ConfidenceLevel.DEFINITIVE)
        assert result.ci_upper - result.ci_lower < 30  # Relatively narrow CI
    
    def test_low_coverage(self):
        """Test with low marker coverage."""
        result = prs_percentile_ci(1.0, n_markers=10, total_markers=100)
        
        assert result.confidence_level in (ConfidenceLevel.LOW, ConfidenceLevel.UNCERTAIN)
        assert "LOW COVERAGE" in str(result.notes)
    
    def test_zero_markers(self):
        """Test with no markers."""
        result = prs_percentile_ci(0, n_markers=0, total_markers=100)
        
        assert result.confidence_level == ConfidenceLevel.UNCERTAIN


class TestAncestrySimilarityStats:
    """Tests for ancestry similarity statistics."""
    
    def test_high_similarity(self):
        """Test high similarity score."""
        result = ancestry_similarity_stats(20, 24, baseline_frequency=0.5)
        
        assert result.value > 80  # High similarity percentage
        assert result.p_value < 0.05  # Significantly different from random
    
    def test_random_similarity(self):
        """Test similarity close to random."""
        result = ancestry_similarity_stats(12, 24, baseline_frequency=0.5)
        
        assert 45 < result.value < 55  # Close to 50%
        assert result.p_value > 0.5  # Not significantly different
    
    def test_zero_markers(self):
        """Test with no markers."""
        result = ancestry_similarity_stats(0, 0, baseline_frequency=0.5)
        
        assert result.confidence_level == ConfidenceLevel.UNCERTAIN


class TestHaplogroupConfidence:
    """Tests for haplogroup confidence calculation."""
    
    def test_definitive_call(self):
        """Test definitive haplogroup call."""
        level, score = haplogroup_confidence(
            defining_snps_found=1,
            defining_snps_total=1,
            consistent_markers=5,
            inconsistent_markers=0
        )
        
        assert level == ConfidenceLevel.DEFINITIVE
        assert score > 0.95
    
    def test_inferred_call(self):
        """Test inferred haplogroup (no defining SNP)."""
        level, score = haplogroup_confidence(
            defining_snps_found=0,
            defining_snps_total=1,
            consistent_markers=4,
            inconsistent_markers=1
        )
        
        assert level in (ConfidenceLevel.MEDIUM, ConfidenceLevel.LOW)
    
    def test_uncertain_call(self):
        """Test uncertain haplogroup."""
        level, score = haplogroup_confidence(
            defining_snps_found=0,
            defining_snps_total=1,
            consistent_markers=1,
            inconsistent_markers=2
        )
        
        assert level in (ConfidenceLevel.LOW, ConfidenceLevel.UNCERTAIN)


class TestDiplotypeConfidence:
    """Tests for diplotype calling confidence."""
    
    def test_high_confidence_diplotype(self):
        """Test high confidence diplotype call."""
        conf, warnings = diplotype_confidence(
            "*4", "*4",
            markers_for_allele_1=3,
            total_markers_allele_1=3,
            markers_for_allele_2=3,
            total_markers_allele_2=3,
            phase_ambiguous=False
        )
        
        assert conf > 0.9
        assert len(warnings) == 0
    
    def test_reference_allele_inference(self):
        """Test when *1 is inferred."""
        conf, warnings = diplotype_confidence(
            "*1", "*4",
            markers_for_allele_1=0,
            total_markers_allele_1=0,  # *1 has no defining markers
            markers_for_allele_2=2,
            total_markers_allele_2=3,
            phase_ambiguous=False
        )
        
        assert any("*1" in w for w in warnings)
        assert conf < 0.95  # Capped due to *1 inference


class TestTraitProbabilityCI:
    """Tests for trait probability confidence interval."""
    
    def test_homozygous_effect(self):
        """Test homozygous effect allele."""
        result = trait_probability_ci(
            genotype="AA",
            effect_allele="A",
            baseline_probability=0.3,
            odds_ratio=2.5,
            odds_ratio_ci_lower=1.8,
            odds_ratio_ci_upper=3.4
        )
        
        # Two effect alleles should increase probability
        assert result.value > 0.3
        assert result.ci_lower < result.value < result.ci_upper
    
    def test_no_effect_allele(self):
        """Test with no effect alleles."""
        result = trait_probability_ci(
            genotype="TT",
            effect_allele="A",
            baseline_probability=0.3,
            odds_ratio=2.5,
            odds_ratio_ci_lower=1.8,
            odds_ratio_ci_upper=3.4
        )
        
        # No effect alleles - should equal baseline
        assert abs(result.value - 0.3) < 0.01


# =============================================================================
# UTILITY FUNCTION TESTS
# =============================================================================

class TestCombineConfidenceLevels:
    """Tests for combining confidence levels."""
    
    def test_takes_minimum(self):
        """Should return minimum confidence level."""
        levels = [ConfidenceLevel.HIGH, ConfidenceLevel.MEDIUM, ConfidenceLevel.LOW]
        result = combine_confidence_levels(levels)
        
        assert result == ConfidenceLevel.LOW
    
    def test_single_level(self):
        """Single level should return itself."""
        result = combine_confidence_levels([ConfidenceLevel.HIGH])
        
        assert result == ConfidenceLevel.HIGH
    
    def test_empty(self):
        """Empty list should return UNCERTAIN."""
        result = combine_confidence_levels([])
        
        assert result == ConfidenceLevel.UNCERTAIN


class TestConfidenceToColor:
    """Tests for confidence level to color mapping."""
    
    def test_all_levels_have_colors(self):
        """All confidence levels should have colors."""
        for level in ConfidenceLevel:
            color = confidence_to_color(level)
            assert color.startswith("#")
            assert len(color) == 7  # #RRGGBB format
    
    def test_high_is_green(self):
        """High confidence should be green-ish."""
        color = confidence_to_color(ConfidenceLevel.HIGH)
        # Green colors typically have high G component
        assert color[3:5] > color[1:3]  # G > R


class TestFormatCIString:
    """Tests for CI string formatting."""
    
    def test_basic_format(self):
        """Test basic formatting."""
        s = format_ci_string(65.0, 58.2, 71.8, decimals=1)
        
        assert "65.0" in s
        assert "58.2" in s
        assert "71.8" in s
    
    def test_percentage_format(self):
        """Test percentage formatting."""
        s = format_ci_string(65.0, 58.2, 71.8, as_percentage=True)
        
        assert "%" in s


# =============================================================================
# INTEGRATION TESTS
# =============================================================================

class TestStatisticalResultIntegration:
    """Integration tests for StatisticalResult."""
    
    def test_result_to_dict(self):
        """Test full result serialization."""
        result = StatisticalResult(
            value=75.0,
            ci_lower=62.0,
            ci_upper=88.0,
            confidence_level=ConfidenceLevel.HIGH,
            p_value=0.001,
            standard_error=8.5,
            n_observations=25,
            n_total=30,
            method="test",
            quality_score=0.83,
            notes=["Test note"]
        )
        
        d = result.to_dict()
        
        assert d["value"] == 75.0
        assert d["confidence"] == "HIGH"
        assert d["p_value"] == 0.001
        assert d["n_markers"] == 25
        assert "notes" in d


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
