"""
Tests for personal-genomics dataset modules.

Run with: pytest tests/test_datasets.py -v
"""

import sys
from pathlib import Path
import pytest

# Add the skill to path
SKILL_PATH = Path(__file__).parent.parent
sys.path.insert(0, str(SKILL_PATH))
sys.path.insert(0, str(SKILL_PATH / "personal_genomics"))


class TestThousandGenomes:
    """Tests for 1000 Genomes dataset."""
    
    def test_import(self):
        """Test that module can be imported."""
        from datasets.thousand_genomes import ThousandGenomes
        assert ThousandGenomes is not None
    
    def test_initialization(self):
        """Test dataset initialization."""
        from datasets.thousand_genomes import ThousandGenomes
        tg = ThousandGenomes()
        assert tg.name == "1000genomes"
        assert tg.version == "phase3_v5b"
    
    def test_download(self):
        """Test dataset download."""
        from datasets.thousand_genomes import ThousandGenomes
        tg = ThousandGenomes()
        result = tg.download()
        assert result is True
        assert tg.is_downloaded
    
    def test_lookup_variant(self):
        """Test variant lookup."""
        from datasets.thousand_genomes import ThousandGenomes
        tg = ThousandGenomes()
        tg.download()
        
        # Look up a known AIM
        var = tg.lookup_variant("rs1426654")
        assert var is not None
        assert var.rsid == "rs1426654"
        assert var.gene == "SLC24A5"
        assert "CEU" in var.frequencies
    
    def test_population_similarity(self):
        """Test population similarity calculation."""
        from datasets.thousand_genomes import ThousandGenomes
        tg = ThousandGenomes()
        tg.download()
        
        # Test with European-like genotypes
        genotypes = {
            "rs1426654": "AA",  # European
            "rs16891982": "GG",  # European
            "rs2814778": "TT",  # Not African
        }
        
        similarities = tg.calculate_population_similarity(genotypes)
        assert len(similarities) > 0
        
        # European populations should score high
        european_pops = ["CEU", "GBR", "FIN", "IBS", "TSI"]
        european_total = sum(similarities.get(p, 0) for p in european_pops)
        assert european_total > 0.5  # Should be >50% European


class TestHGDP:
    """Tests for HGDP dataset."""
    
    def test_import(self):
        """Test that module can be imported."""
        from datasets.hgdp import HGDP
        assert HGDP is not None
    
    def test_initialization(self):
        """Test dataset initialization."""
        from datasets.hgdp import HGDP
        hgdp = HGDP()
        assert hgdp.name == "hgdp"
    
    def test_populations_defined(self):
        """Test that population metadata is defined."""
        from datasets.hgdp import HGDP_POPULATIONS, HGDP_REGIONS
        assert len(HGDP_POPULATIONS) > 0
        assert len(HGDP_REGIONS) > 0
        assert "Europe" in HGDP_REGIONS


class TestSGDP:
    """Tests for SGDP dataset."""
    
    def test_import(self):
        """Test that module can be imported."""
        from datasets.sgdp import SGDPDataset
        assert SGDPDataset is not None
    
    def test_populations_defined(self):
        """Test that population metadata is defined."""
        from datasets.sgdp import SGDP_POPULATIONS, ALL_SGDP_POPULATIONS
        assert len(ALL_SGDP_POPULATIONS) > 100
        assert "Europe" in SGDP_POPULATIONS


class TestAncientDNA:
    """Tests for Ancient DNA dataset."""
    
    def test_import(self):
        """Test that module can be imported."""
        from datasets.ancient_dna import AncientDNADataset
        assert AncientDNADataset is not None
    
    def test_markers_defined(self):
        """Test that ancient markers are defined."""
        from datasets.ancient_dna import ANCIENT_MARKERS, MARKERS_BY_POPULATION
        assert len(ANCIENT_MARKERS) > 0
        assert "WHG" in MARKERS_BY_POPULATION
        assert "Neanderthal" in MARKERS_BY_POPULATION
    
    def test_ancestral_signal_analysis(self):
        """Test ancestral signal analysis."""
        from datasets.ancient_dna import AncientDNADataset
        
        ds = AncientDNADataset()
        ds.download()
        
        # Test with European-like genotypes
        genotypes = {
            "rs12913832": "GG",  # Blue eyes
            "rs1426654": "AA",  # Light skin
        }
        
        results = ds.analyze_ancestral_signals(genotypes)
        assert "WHG" in results
        assert "Neanderthal" in results


class TestClinVar:
    """Tests for ClinVar dataset."""
    
    def test_import(self):
        """Test that module can be imported."""
        from datasets.clinvar import ClinVar
        assert ClinVar is not None
    
    def test_download(self):
        """Test dataset download."""
        from datasets.clinvar import ClinVar
        cv = ClinVar()
        result = cv.download()
        assert result is True
    
    def test_pathogenic_lookup(self):
        """Test looking up pathogenic variants."""
        from datasets.clinvar import ClinVar
        cv = ClinVar()
        cv.download()
        
        # Look up a known pathogenic variant
        var = cv.get_clinvar_annotation("rs6025")
        assert var is not None
        assert "pathogenic" in var.clinical_significance.lower()


class TestPharmGKB:
    """Tests for PharmGKB dataset."""
    
    def test_import(self):
        """Test that module can be imported."""
        from datasets.pharmgkb import PharmGKB
        assert PharmGKB is not None
    
    def test_cpic_guidelines(self):
        """Test CPIC guidelines loaded."""
        from datasets.pharmgkb import PharmGKB, CPIC_GUIDELINES
        assert len(CPIC_GUIDELINES) > 0
        assert "CYP2D6" in CPIC_GUIDELINES
    
    def test_drug_interactions(self):
        """Test getting drug interactions."""
        from datasets.pharmgkb import PharmGKB
        pgkb = PharmGKB()
        pgkb.download()
        
        interactions = pgkb.get_drug_interactions("CYP2D6")
        assert len(interactions) > 0


class TestPGSCatalog:
    """Tests for PGS Catalog dataset."""
    
    def test_import(self):
        """Test that module can be imported."""
        from datasets.pgs_catalog import PGSCatalog
        assert PGSCatalog is not None
    
    def test_available_prs(self):
        """Test getting available PRS models."""
        from datasets.pgs_catalog import PGSCatalog
        pgs = PGSCatalog()
        pgs.download()
        
        models = pgs.get_available_prs()
        assert len(models) > 0
        
        # Check that key conditions are covered
        traits = [m.trait for m in models]
        assert any("coronary" in t.lower() for t in traits)
        assert any("diabetes" in t.lower() for t in traits)


class TestGWASCatalog:
    """Tests for GWAS Catalog dataset."""
    
    def test_import(self):
        """Test that module can be imported."""
        from datasets.gwas_catalog import GWASCatalog
        assert GWASCatalog is not None
    
    def test_download(self):
        """Test dataset download."""
        from datasets.gwas_catalog import GWASCatalog
        gwas = GWASCatalog()
        result = gwas.download()
        assert result is True
    
    def test_trait_lookup(self):
        """Test looking up traits."""
        from datasets.gwas_catalog import GWASCatalog
        gwas = GWASCatalog()
        gwas.download()
        
        # Look up a known variant-trait association
        associations = gwas.get_variant_associations("rs12913832")
        assert len(associations) > 0
        assert any("eye" in a.trait.lower() for a in associations)


class TestGnomAD:
    """Tests for gnomAD dataset."""
    
    def test_import(self):
        """Test that module can be imported."""
        from datasets.gnomad import GnomAD
        assert GnomAD is not None
    
    def test_download(self):
        """Test dataset download."""
        from datasets.gnomad import GnomAD
        gnomad = GnomAD()
        result = gnomad.download()
        assert result is True
    
    def test_frequency_lookup(self):
        """Test allele frequency lookup."""
        from datasets.gnomad import GnomAD
        gnomad = GnomAD()
        gnomad.download()
        
        # Look up a known common variant
        freq = gnomad.get_allele_frequency("rs429358")
        assert freq is not None
        assert 0 < freq < 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
