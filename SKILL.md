# Personal Genomics Skill v3.0

Comprehensive local DNA analysis with **1300+ markers** across **16 categories**. Privacy-first genetic analysis for AI agents.

## Quick Start

```bash
python comprehensive_analysis.py /path/to/dna_file.txt
```

## Triggers

Activate this skill when user mentions:
- DNA analysis, genetic analysis, genome analysis
- 23andMe, AncestryDNA, MyHeritage results
- Pharmacogenomics, drug-gene interactions
- Genetic risk, disease risk, health risk
- Carrier status, carrier testing
- VCF file analysis
- APOE, MTHFR, CYP2D6, or other gene names
- Polygenic risk scores

## Supported Files

- 23andMe, AncestryDNA, MyHeritage, FTDNA
- VCF files (whole genome/exome, .vcf or .vcf.gz)
- Any tab-delimited rsid format

## Output Location

`~/dna-analysis/reports/`

- `agent_summary.json` - AI-optimized, priority-sorted
- `full_analysis.json` - Complete data
- `report.txt` - Human-readable

## Marker Categories (16 total)

1. **Pharmacogenomics** (159) - Drug metabolism (CYP450, DPYD, TPMT, HLA)
2. **Polygenic Risk Scores** (277) - CAD, T2D, cancer, Alzheimer's, etc.
3. **Carrier Status** (140) - CF, sickle cell, Tay-Sachs, 35+ diseases
4. **Health Risks** (213) - APOE, Factor V Leiden, AMD, celiac
5. **Traits** (163) - Eye color, taste, alcohol flush, chronotype
6. **Nutrition** (34) - MTHFR, vitamin D, caffeine, omega-3
7. **Fitness** (30) - ACTN3, injury risk, recovery
8. **Neurogenetics** (28) - COMT, BDNF, dopamine, serotonin
9. **Longevity** (30) - FOXO3, TERT, aging markers
10. **Immunity** (43) - HLA, autoimmunity, infection risk
11. **Rare Diseases** (29) - Lysosomal storage, connective tissue, neurological
12. **Mental Health** (25) - Depression, anxiety, ADHD, substance use
13. **Dermatology** (37) - Melanoma, psoriasis, eczema, pigmentation
14. **Vision & Hearing** (33) - AMD, glaucoma, hearing loss, ototoxicity
15. **Fertility** (31) - PCOS, endometriosis, pregnancy risk
16. **Ancestry** (78) - Population informative markers

## Agent Integration

The `agent_summary.json` provides:

```json
{
  "critical_alerts": [],      // MUST share with healthcare
  "high_priority": [],        // Important findings
  "medium_priority": [],      // Moderate relevance
  "pharmacogenomics_alerts": [],
  "apoe_status": {},
  "polygenic_risk_scores": {},
  "lifestyle_recommendations": {
    "diet": [],
    "exercise": [],
    "supplements": [],
    "avoid": []
  },
  "drug_interaction_matrix": {},
  "notable_traits": []
}
```

## Critical Findings (Always Alert User)

### Pharmacogenomics
- **DPYD** variants - 5-FU/capecitabine FATAL toxicity risk
- **HLA-B*5701** - Abacavir hypersensitivity
- **HLA-B*1502** - Carbamazepine SJS (certain populations)
- **MT-RNR1** - Aminoglycoside-induced deafness

### Disease Risk
- **APOE ε4/ε4** - ~12x Alzheimer's risk
- **Factor V Leiden** - Thrombosis risk, contraceptive implications
- **BRCA proxies** - Breast cancer risk

### Carrier Status
- **CFTR** - Cystic fibrosis (1 in 25 Europeans)
- **HBB** - Sickle cell (1 in 12 African Americans)
- **HEXA** - Tay-Sachs (1 in 30 Ashkenazi Jews)

## Interpretation Resources

See `references/INTERPRETATION_GUIDE.md` for:
- How to explain findings in layperson language
- Risk communication best practices
- When to recommend genetic counseling
- Cultural sensitivity guidelines

## Privacy

- All analysis runs locally
- Zero network requests
- No data leaves the machine

## Limitations

- Consumer arrays miss rare variants (~0.1% of genome)
- Results are probabilistic, not deterministic
- Not a medical diagnosis
- Most conditions 50-80% non-genetic
- Consult healthcare providers for medical decisions

## When to Recommend Genetic Counseling

- APOE ε4/ε4 genotype
- Multiple critical pharmacogenomic findings
- Carrier status with reproduction implications
- Hereditary cancer syndrome indicators
- Results causing significant user distress
