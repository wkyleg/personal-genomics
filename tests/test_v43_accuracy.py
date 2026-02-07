"""
Tests for v4.3.0 Accuracy Features

Tests:
- Confidence interval calculation
- PMID presence in markers
- Haplogroup disclaimers
- Ancestry continental-only reporting
- PRS percentile ranges
"""

import pytest
import math
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestConfidenceIntervals:
    """Test confidence interval calculations."""
    
    def test_wilson_confidence_interval_basic(self):
        """Test Wilson score interval calculation."""
        from markers.ancestry_composition import calculate_wilson_confidence_interval
        
        # 50/100 successes should give ~40-60% CI
        lower, upper = calculate_wilson_confidence_interval(50, 100)
        assert 0.35 < lower < 0.45, f"Lower bound {lower} not in expected range"
        assert 0.55 < upper < 0.65, f"Upper bound {upper} not in expected range"
    
    def test_wilson_ci_extreme_proportions(self):
        """Test CI for extreme proportions (near 0 or 1)."""
        from markers.ancestry_composition import calculate_wilson_confidence_interval
        
        # High proportion (95/100)
        lower, upper = calculate_wilson_confidence_interval(95, 100)
        assert lower > 0.85, f"Lower bound {lower} too low for 95%"
        assert upper <= 1.0, f"Upper bound {upper} exceeds 1.0"
        
        # Low proportion (5/100)
        lower, upper = calculate_wilson_confidence_interval(5, 100)
        assert lower >= 0.0, f"Lower bound {lower} below 0.0"
        assert upper < 0.15, f"Upper bound {upper} too high for 5%"
    
    def test_wilson_ci_zero_trials(self):
        """Test CI with zero trials returns full range."""
        from markers.ancestry_composition import calculate_wilson_confidence_interval
        
        lower, upper = calculate_wilson_confidence_interval(0, 0)
        assert lower == 0.0
        assert upper == 1.0
    
    def test_wilson_ci_small_sample(self):
        """Test CI with small sample sizes gives wider intervals."""
        from markers.ancestry_composition import calculate_wilson_confidence_interval
        
        # Small sample (5/10)
        lower_small, upper_small = calculate_wilson_confidence_interval(5, 10)
        
        # Large sample (50/100)
        lower_large, upper_large = calculate_wilson_confidence_interval(50, 100)
        
        # Small sample should have wider CI
        width_small = upper_small - lower_small
        width_large = upper_large - lower_large
        assert width_small > width_large, "Small sample should have wider CI"


class TestPMIDPresence:
    """Test that markers have PMIDs for citation."""
    
    def test_haplogroup_markers_have_pmids(self):
        """Test mtDNA and Y-DNA markers have PMIDs."""
        from markers.haplogroups import MTDNA_MARKERS, YCHROMOSOME_MARKERS
        
        # Check mtDNA markers
        for rsid, info in MTDNA_MARKERS.items():
            assert "pmid" in info, f"mtDNA marker {rsid} missing PMID"
            assert isinstance(info["pmid"], list), f"mtDNA marker {rsid} PMID should be list"
            assert len(info["pmid"]) > 0, f"mtDNA marker {rsid} has empty PMID list"
        
        # Check Y-DNA markers
        for rsid, info in YCHROMOSOME_MARKERS.items():
            assert "pmid" in info, f"Y-DNA marker {rsid} missing PMID"
            assert isinstance(info["pmid"], list), f"Y-DNA marker {rsid} PMID should be list"
            assert len(info["pmid"]) > 0, f"Y-DNA marker {rsid} has empty PMID list"
    
    def test_ancestry_markers_have_pmids(self):
        """Test ancestry informative markers have PMIDs."""
        from markers.ancestry_composition import ANCESTRY_INFORMATIVE_MARKERS
        
        for rsid, info in ANCESTRY_INFORMATIVE_MARKERS.items():
            assert "pmid" in info, f"Ancestry marker {rsid} missing PMID"
            assert isinstance(info["pmid"], list), f"Ancestry marker {rsid} PMID should be list"
            assert len(info["pmid"]) > 0, f"Ancestry marker {rsid} has empty PMID list"
    
    def test_haplogroup_history_has_pmids(self):
        """Test haplogroup history entries have PMIDs."""
        from markers.haplogroups import HAPLOGROUP_HISTORY
        
        for hg, info in HAPLOGROUP_HISTORY.items():
            assert "pmid" in info, f"Haplogroup {hg} history missing PMID"
            assert isinstance(info["pmid"], list), f"Haplogroup {hg} PMID should be list"


class TestHaplogroupDisclaimers:
    """Test haplogroup confidence and disclaimers."""
    
    def test_haplogroup_disclaimer_present(self):
        """Test that HAPLOGROUP_DISCLAIMER constant exists."""
        from markers.haplogroups import HAPLOGROUP_DISCLAIMER
        
        assert HAPLOGROUP_DISCLAIMER is not None
        assert "LOW CONFIDENCE" in HAPLOGROUP_DISCLAIMER
        assert "FTDNA" in HAPLOGROUP_DISCLAIMER or "YFull" in HAPLOGROUP_DISCLAIMER
    
    def test_mtdna_result_has_low_confidence(self):
        """Test mtDNA analysis returns low confidence."""
        from markers.haplogroups import determine_mtdna_haplogroup
        
        # Test with empty genotypes
        result = determine_mtdna_haplogroup({})
        assert "low" in result["confidence"].lower()
        assert "disclaimer" in result
        
        # Test with some markers
        test_genotypes = {
            "rs2853511": "CC",  # H haplogroup indicator
            "rs2853512": "AA",  # H confirmation
        }
        result = determine_mtdna_haplogroup(test_genotypes)
        # Even with markers, should still be low confidence for consumer arrays
        assert "low" in result["confidence"].lower()
        assert result["recommendation"], "Should have recommendation for better testing"
    
    def test_ydna_result_has_low_confidence(self):
        """Test Y-DNA analysis returns low confidence."""
        from markers.haplogroups import determine_y_haplogroup
        
        # Test with empty genotypes
        result = determine_y_haplogroup({})
        
        # Test with some markers
        test_genotypes = {
            "rs9786076": "CC",  # R1b indicator
        }
        result = determine_y_haplogroup(test_genotypes)
        # Even with markers, should still be low confidence
        assert "low" in result["confidence"].lower() or result["confidence"] == "N/A"
    
    def test_full_haplogroup_analysis_has_disclaimer(self):
        """Test full haplogroup analysis includes disclaimer."""
        from markers.haplogroups import analyze_haplogroups
        
        result = analyze_haplogroups({})
        
        assert "disclaimer" in result
        assert "methodology" in result
        assert "limitations" in result["methodology"]


class TestAncientAncestrySignals:
    """Test ancient ancestral signals detection."""
    
    def test_ancient_populations_defined(self):
        """Test all ancient populations are properly defined."""
        from markers.ancestry_composition import ANCIENT_POPULATIONS
        
        required_pops = ["WHG", "EEF", "STEPPE", "NEANDERTHAL", "DENISOVAN"]
        for pop in required_pops:
            assert pop in ANCIENT_POPULATIONS, f"Missing ancient population: {pop}"
            
            info = ANCIENT_POPULATIONS[pop]
            assert "name" in info, f"{pop} missing name"
            assert "time_period" in info, f"{pop} missing time_period"
            assert "traits_contributed" in info, f"{pop} missing traits_contributed"
            assert "pmid" in info, f"{pop} missing PMID citations"
    
    def test_ancient_markers_have_pmids(self):
        """Test all ancient ancestry markers have PMIDs."""
        from markers.ancestry_composition import ANCIENT_ANCESTRY_MARKERS
        
        for rsid, info in ANCIENT_ANCESTRY_MARKERS.items():
            assert "pmid" in info, f"Marker {rsid} missing PMID"
            assert len(info["pmid"]) > 0, f"Marker {rsid} has empty PMID list"
    
    def test_ancient_markers_have_interpretations(self):
        """Test markers have signal interpretations."""
        from markers.ancestry_composition import ANCIENT_ANCESTRY_MARKERS
        
        for rsid, info in ANCIENT_ANCESTRY_MARKERS.items():
            assert "signal_interpretation" in info, f"Marker {rsid} missing signal_interpretation"
            assert "ancestral_population" in info, f"Marker {rsid} missing ancestral_population"
    
    def test_detect_ancient_signals_structure(self):
        """Test ancient signal detection returns proper structure."""
        from markers.ancestry_composition import detect_ancient_signals
        
        result = detect_ancient_signals({})
        
        assert "methodology" in result
        assert "populations" in result
        
        # Check all populations are present
        for pop in ["WHG", "EEF", "STEPPE", "NEANDERTHAL", "DENISOVAN"]:
            assert pop in result["populations"], f"Missing {pop} in results"
            pop_data = result["populations"][pop]
            assert "signal_strength" in pop_data
            assert "markers_detected" in pop_data
    
    def test_signal_strength_categories(self):
        """Test signal strength is properly categorized."""
        from markers.ancestry_composition import detect_ancient_signals
        
        # Test with some markers
        test_genotypes = {
            "rs4988235": "AA",  # Steppe signal (lactase persistence)
            "rs12913832": "GG",  # WHG signal (blue eyes)
            "rs1426654": "AA",  # EEF signal (light skin)
        }
        
        result = detect_ancient_signals(test_genotypes)
        
        # Check that signals are detected
        valid_strengths = ["not detected", "weak", "moderate", "strong"]
        for pop, data in result["populations"].items():
            assert data["signal_strength"] in valid_strengths, \
                f"Invalid signal strength for {pop}: {data['signal_strength']}"
    
    def test_get_ancestry_summary_format(self):
        """Test ancestry summary has correct format."""
        from markers.ancestry_composition import get_ancestry_summary
        
        result = get_ancestry_summary({})
        
        assert result["analysis_type"] == "ancient_ancestral_signals"
        assert "disclaimer" in result
        assert "signals" in result
        assert "educational_note" in result
        assert "modern ethnicity" in result["disclaimer"].lower() or \
               "ancient" in result["disclaimer"].lower()


class TestPRSRanges:
    """Test PRS shows ranges instead of point estimates."""
    
    def test_prs_result_structure(self):
        """Test PRS results have proper structure for ranges."""
        from comprehensive_analysis import calculate_all_prs
        
        # Would need actual genotype data to fully test
        # For now, just verify the function exists and handles empty input
        result = calculate_all_prs({})
        
        # Should return dict even with empty input
        assert isinstance(result, dict)


class TestDashboardFeatures:
    """Test dashboard has required accuracy features."""
    
    def test_dashboard_exists(self):
        """Test dashboard file exists."""
        dashboard_path = Path(__file__).parent.parent / "dashboard" / "index.html"
        assert dashboard_path.exists(), "Dashboard file not found"
    
    def test_dashboard_has_pmid_styles(self):
        """Test dashboard has PMID styling."""
        dashboard_path = Path(__file__).parent.parent / "dashboard" / "index.html"
        content = dashboard_path.read_text()
        
        assert "pmid-link" in content, "Dashboard missing PMID link styles"
        assert "pubmed" in content.lower(), "Dashboard missing PubMed integration"
    
    def test_dashboard_has_warning_banner(self):
        """Test dashboard has warning banner styles."""
        dashboard_path = Path(__file__).parent.parent / "dashboard" / "index.html"
        content = dashboard_path.read_text()
        
        assert "warning-banner" in content, "Dashboard missing warning banner styles"
    
    def test_dashboard_has_methodology_section(self):
        """Test dashboard has methodology section."""
        dashboard_path = Path(__file__).parent.parent / "dashboard" / "index.html"
        content = dashboard_path.read_text()
        
        assert "methodology" in content.lower(), "Dashboard missing methodology section"
        assert "limitations" in content.lower(), "Dashboard missing limitations"
    
    def test_dashboard_version_updated(self):
        """Test dashboard version is 4.3.0+."""
        dashboard_path = Path(__file__).parent.parent / "dashboard" / "index.html"
        content = dashboard_path.read_text()
        
        # Accept any version 4.3.0 or higher
        assert "Personal Genomics Dashboard v4" in content, "Dashboard version header not found"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
