# Personal Genomics v3.0

Privacy-first local DNA analysis for AI agents. Comprehensive genetic analysis with **1300+ validated markers** across **16 categories** covering pharmacogenomics, disease risk, carrier status, traits, mental health, rare diseases, and lifestyle factors.

![Logo](logo.svg)

## Features

### Core Analysis
- **1300+ validated genetic markers** across 16 categories
- **Polygenic Risk Scores (PRS)** for 10+ major conditions with confidence intervals
- **Pharmacogenomics** with CPIC Level 1A drug-gene interactions
- **Carrier screening** for 35+ recessive diseases including rare diseases
- **VCF support** for whole genome/exome sequencing
- **Agent-friendly JSON output** with priority-sorted actionable items
- **Zero network requests** - all analysis runs locally
- **Universal ethnic coverage** - works for all ancestries worldwide

### New in v3.0
- 🧬 **Rare Diseases** - Lysosomal storage disorders, connective tissue, neurological conditions
- 🧠 **Mental Health** - Depression, anxiety, bipolar, ADHD, substance use genetics
- 🌞 **Dermatology** - Skin cancer risk, psoriasis, eczema, pigmentation, aging
- 👁️ **Vision & Hearing** - AMD, glaucoma, hearing loss, ototoxicity
- 🤰 **Fertility** - PCOS, endometriosis, male fertility, pregnancy complications
- 💊 **Drug Interaction Matrix** - Critical interactions, warnings, dosing adjustments
- 🏃 **Lifestyle Recommendations Engine** - Personalized diet, exercise, supplement suggestions
- 📚 **Interpretation Guide** - How to communicate genetic findings responsibly

## Supported Formats

- 23andMe (v3, v4, v5)
- AncestryDNA
- MyHeritage
- FamilyTreeDNA
- Nebula Genomics
- VCF files (whole genome/exome, gzipped supported)
- Any tab-delimited rsid format

## Installation

```bash
# Install via clawhub (recommended)
clawhub install personal-genomics

# Or clone directly
git clone https://github.com/wkyleg/personal-genomics.git
```

## Usage

### Command Line

```bash
python comprehensive_analysis.py /path/to/dna_file.txt
```

### As OpenClaw Skill

```
Analyze my DNA file at ~/Downloads/genome.txt
```

### Output Files

Reports are saved to `~/dna-analysis/reports/`:

- `agent_summary.json` - AI-optimized output with priority-sorted actionable items
- `full_analysis.json` - Complete analysis data
- `report.txt` - Human-readable report

## Marker Categories

| Category | Markers | Description |
|----------|---------|-------------|
| Pharmacogenomics | 159 | Drug metabolism (CYP450, DPYD, TPMT, HLA) |
| Polygenic Risk | 277 | Disease risk scores (CAD, T2D, cancer, etc.) |
| Carrier Status | 140 | Recessive disease carriers (CF, sickle cell, Tay-Sachs) |
| Health Risks | 213 | Disease susceptibility (APOE, Factor V, AMD) |
| Traits | 163 | Physical, sensory, behavioral traits |
| Rare Diseases | 29 | Rare genetic conditions (lysosomal, neurological) |
| Mental Health | 25 | Psychiatric genetics (BDNF, COMT, FKBP5) |
| Dermatology | 37 | Skin conditions (MC1R, FLG, psoriasis) |
| Vision & Hearing | 33 | Sensory conditions (AMD, glaucoma, hearing loss) |
| Fertility | 31 | Reproductive health (PCOS, pregnancy risk) |
| Nutrition | 34 | Nutrigenomics (MTHFR, vitamin D, caffeine) |
| Fitness | 30 | Athletic performance, injury risk, recovery |
| Neurogenetics | 28 | Cognition, behavior, mental health |
| Longevity | 30 | Aging and lifespan markers |
| Immunity | 43 | HLA, autoimmunity, infection susceptibility |
| Ancestry | 78 | Population informative markers |

## Agent-Friendly Output

The `agent_summary.json` is designed for AI agents to quickly identify what matters:

```json
{
  "critical_alerts": [...],
  "high_priority": [...],
  "medium_priority": [...],
  "pharmacogenomics_alerts": [...],
  "apoe_status": {
    "genotype": "ε3/ε4",
    "risk_level": "elevated",
    "recommendations": [...]
  },
  "polygenic_risk_scores": {
    "cad": {"percentile_estimate": 75, "confidence": "moderate", "confidence_interval": [65, 85]},
    "t2d": {"percentile_estimate": 42, "confidence": "moderate"}
  },
  "lifestyle_recommendations": {
    "diet": [...],
    "exercise": [...],
    "supplements": [...],
    "avoid": [...]
  },
  "drug_interaction_matrix": {
    "critical_interactions": [...],
    "warnings": [...],
    "dosing_adjustments": [...]
  }
}
```

## Critical Pharmacogenomics

| Gene | Drugs Affected | Clinical Impact |
|------|----------------|-----------------|
| DPYD | 5-FU, capecitabine | **Fatal toxicity risk** |
| HLA-B*5701 | Abacavir | Hypersensitivity |
| HLA-B*1502 | Carbamazepine | Stevens-Johnson Syndrome |
| MT-RNR1 | Aminoglycosides | Irreversible deafness |
| CYP2D6 | Codeine, tramadol | Prodrug activation |
| CYP2C19 | Clopidogrel (Plavix) | Platelet inhibition |
| SLCO1B1 | Simvastatin | Myopathy risk |

## Polygenic Risk Score Conditions

- Coronary Artery Disease (CAD)
- Type 2 Diabetes (T2D)
- Breast Cancer
- Prostate Cancer
- Colorectal Cancer
- Alzheimer's Disease
- Atrial Fibrillation
- Inflammatory Bowel Disease
- Obesity
- Major Depression

## Interpretation Guide

See `references/INTERPRETATION_GUIDE.md` for:
- How to explain findings to users in layperson language
- Risk communication best practices
- When to recommend genetic counseling
- Cultural sensitivity in ancestry discussions
- Handling emotional responses to results

## Example Outputs

See `examples/` directory for:
- `agent_summary.json` - Typical analysis output
- `sample_high_risk.json` - High-risk individual (critical alerts)
- `sample_carrier.json` - Multiple carrier status findings

## Testing

```bash
# Install test dependencies
pip install pytest

# Run all tests
pytest tests/ -v
```

74 comprehensive tests covering all marker modules, VCF parsing, and edge cases.

## Privacy

- **All analysis runs locally** - no data leaves your machine
- **No network requests** - completely offline capable
- **No tracking or telemetry**
- **Your genetic data stays yours**

## Limitations

1. **Not diagnostic** - Results are informational, not medical diagnoses
2. **Array limitations** - Consumer arrays capture ~0.1% of genome; rare variants often missed
3. **Probabilistic** - Polygenic scores indicate risk, not certainty
4. **Environment matters** - Most conditions are 50-80% non-genetic
5. **Population effects** - Some markers better validated in European ancestry
6. **No structural variants** - CNVs and large deletions not detected

## Contributing

Contributions welcome! Please ensure new markers include:

- rsID
- Gene name
- Risk/effect allele
- Evidence level (strong/moderate/preliminary)
- Source citation (PMID or ClinVar ID)
- Actionable recommendations (if applicable)

## Data Sources

- **PharmGKB** - Pharmacogenomics annotations
- **ClinVar** - Clinical variant interpretations
- **NHGRI-EBI GWAS Catalog** - Genome-wide associations
- **CPIC** - Clinical Pharmacogenetics Implementation Consortium
- **PGS Catalog** - Polygenic Score database
- **OMIM** - Rare disease genetics

## References

1. Relling MV, Klein TE. CPIC: Clinical Pharmacogenetics Implementation Consortium. Clin Pharmacol Ther. 2011.
2. Landrum MJ, et al. ClinVar: public archive of clinically relevant variants. Nucleic Acids Res. 2016.
3. Buniello A, et al. The NHGRI-EBI GWAS Catalog. Nucleic Acids Res. 2019.
4. Whirl-Carrillo M, et al. Pharmacogenomics knowledge for personalized medicine. Clin Pharmacol Ther. 2012.
5. Lambert SA, et al. The Polygenic Score Catalog. Nat Genet. 2021.

## License

MIT License - See LICENSE file for details.

---

**Disclaimer**: This tool is for educational and research purposes. It is not a substitute for professional medical advice, diagnosis, or treatment. Always consult qualified healthcare providers for medical decisions. Critical pharmacogenomic findings should be confirmed with clinical-grade testing before making treatment decisions.
