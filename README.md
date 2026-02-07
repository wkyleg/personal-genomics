# Personal Genomics

Privacy-first local DNA analysis for AI agents. Comprehensive genetic analysis with 2000+ markers covering pharmacogenomics, disease risk, carrier status, traits, and lifestyle factors.

![Logo](logo.svg)

## Features

- **2000+ validated genetic markers** across 10 categories
- **Polygenic Risk Scores** for 10 major conditions (CAD, T2D, cancer, Alzheimer's, etc.)
- **Pharmacogenomics** with CPIC Level 1A drug-gene interactions
- **Carrier screening** for 35+ recessive diseases
- **VCF support** for whole genome/exome sequencing
- **Agent-friendly JSON output** with priority-sorted actionable items
- **Zero network requests** - all analysis runs locally
- **Universal ethnic coverage** - works for all ancestries worldwide

## Supported Formats

- 23andMe (v3, v4, v5)
- AncestryDNA
- MyHeritage
- FamilyTreeDNA
- Nebula Genomics
- VCF files (whole genome/exome)
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
| Pharmacogenomics | 100+ | Drug metabolism (CYP450, DPYD, TPMT, HLA) |
| Polygenic Risk | 150+ | Disease risk scores (CAD, T2D, cancer, etc.) |
| Carrier Status | 35+ | Recessive disease carriers (CF, sickle cell, Tay-Sachs) |
| Health Risks | 50+ | Disease susceptibility (APOE, Factor V, AMD) |
| Traits | 60+ | Physical, sensory, behavioral traits |
| Nutrition | 40+ | Nutrigenomics (MTHFR, vitamin D, caffeine) |
| Fitness | 35+ | Athletic performance, injury risk, recovery |
| Neurogenetics | 35+ | Cognition, behavior, mental health |
| Longevity | 30+ | Aging and lifespan markers |
| Immunity | 50+ | HLA, autoimmunity, infection susceptibility |

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
    "cad": {"percentile_estimate": 75, "confidence": "moderate"},
    "t2d": {"percentile_estimate": 42, "confidence": "moderate"}
  }
}
```

## Pharmacogenomics Coverage

Critical drug-gene pairs from CPIC guidelines:

| Gene | Drugs Affected | Clinical Impact |
|------|----------------|-----------------|
| CYP2D6 | Codeine, tramadol, tamoxifen | Prodrug activation |
| CYP2C19 | Clopidogrel (Plavix), PPIs | Platelet inhibition |
| CYP2C9 + VKORC1 | Warfarin | Bleeding risk |
| DPYD | 5-FU, capecitabine | **Fatal toxicity risk** |
| TPMT/NUDT15 | Azathioprine, 6-MP | Myelosuppression |
| HLA-B*5701 | Abacavir | Hypersensitivity |
| HLA-B*1502 | Carbamazepine | Stevens-Johnson Syndrome |
| SLCO1B1 | Simvastatin | Myopathy risk |

## Polygenic Risk Score Conditions

Validated PRS models for:

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

## Data Sources

Markers are sourced from peer-reviewed databases:

- **PharmGKB** - Pharmacogenomics annotations
- **ClinVar** - Clinical variant interpretations
- **NHGRI-EBI GWAS Catalog** - Genome-wide associations
- **CPIC** - Clinical Pharmacogenetics Implementation Consortium
- **PGS Catalog** - Polygenic Score database

## Privacy

- **All analysis runs locally** - no data leaves your machine
- **No network requests** - completely offline capable
- **No tracking or telemetry**
- **Your genetic data stays yours**

## Limitations

1. **Not diagnostic** - Results are informational, not medical diagnoses
2. **Array limitations** - Consumer arrays capture ~0.1% of genome; rare variants missed
3. **Probabilistic** - Polygenic scores indicate risk, not certainty
4. **Environment matters** - Most conditions are 50-80% non-genetic
5. **Population effects** - Some markers better validated in European ancestry

## Contributing

Contributions welcome! Please ensure new markers include:

- rsID
- Gene name
- Risk/effect allele
- Evidence level
- Source citation
- Actionable recommendations (if applicable)

## References

1. Purcell S, et al. PLINK: a tool set for whole-genome association. Am J Hum Genet. 2007;81(3):559-575.
2. Landrum MJ, et al. ClinVar: public archive of interpretations of clinically relevant variants. Nucleic Acids Res. 2016;44(D1):D862-8.
3. Whirl-Carrillo M, et al. Pharmacogenomics knowledge for personalized medicine. Clin Pharmacol Ther. 2012;92(4):414-7.
4. Relling MV, Klein TE. CPIC: Clinical Pharmacogenetics Implementation Consortium. Clin Pharmacol Ther. 2011;89(3):464-7.
5. Buniello A, et al. The NHGRI-EBI GWAS Catalog. Nucleic Acids Res. 2019;47(D1):D1005-D1012.

## License

MIT License - See LICENSE file for details.

---

**Disclaimer**: This tool is for educational and research purposes. It is not a substitute for professional medical advice, diagnosis, or treatment. Always consult qualified healthcare providers for medical decisions.
