# Personal Genomics Skill v4.0

Comprehensive local DNA analysis with **1450+ markers** across **21 categories**. Privacy-first genetic analysis for AI agents.

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
- APOE, MTHFR, CYP2D6, BRCA, or other gene names
- Polygenic risk scores
- Haplogroups, maternal lineage, paternal lineage
- Ancestry composition, ethnicity
- Hereditary cancer, Lynch syndrome
- Autoimmune genetics, HLA, celiac
- Pain sensitivity, opioid response

## Supported Files

- 23andMe, AncestryDNA, MyHeritage, FTDNA
- VCF files (whole genome/exome, .vcf or .vcf.gz)
- Any tab-delimited rsid format

## Output Location

`~/dna-analysis/reports/`

- `agent_summary.json` - AI-optimized, priority-sorted
- `full_analysis.json` - Complete data
- `report.txt` - Human-readable
- `genetic_report.pdf` - Professional PDF report

## New v4.0 Features

### Haplogroup Analysis
- Mitochondrial DNA (mtDNA) - maternal lineage
- Y-chromosome - paternal lineage (males only)
- Migration history context
- PhyloTree/ISOGG standards

### Ancestry Composition
- Population comparisons (EUR, AFR, EAS, SAS, AMR)
- Admixture detection
- Ancestry informative markers

### Hereditary Cancer Panel
- BRCA1/BRCA2 comprehensive
- Lynch syndrome (MLH1, MSH2, MSH6, PMS2)
- Other genes (APC, TP53, CHEK2, PALB2, ATM)
- ACMG-style classification

### Autoimmune HLA
- Celiac (DQ2/DQ8) - can rule out if negative
- Type 1 Diabetes
- Ankylosing spondylitis (HLA-B27)
- Rheumatoid arthritis, lupus, MS

### Pain Sensitivity
- COMT Val158Met
- OPRM1 opioid receptor
- SCN9A pain signaling
- TRPV1 capsaicin sensitivity
- Migraine susceptibility

### PDF Reports
- Professional format
- Physician-shareable
- Executive summary
- Detailed findings
- Disclaimers included

### Data Quality
- Call rate analysis
- Platform detection
- Confidence scoring
- Quality warnings

### Export Formats
- Genetic counselor clinical export
- Apple Health compatible
- API-ready JSON
- Integration hooks

## Marker Categories (21 total)

1. **Pharmacogenomics** (159) - Drug metabolism
2. **Polygenic Risk Scores** (277) - Disease risk
3. **Carrier Status** (181) - Recessive carriers
4. **Health Risks** (233) - Disease susceptibility
5. **Traits** (163) - Physical/behavioral
6. **Haplogroups** (44) - Lineage markers
7. **Ancestry** (124) - Population informative
8. **Hereditary Cancer** (41) - BRCA, Lynch, etc.
9. **Autoimmune HLA** (31) - HLA associations
10. **Pain Sensitivity** (20) - Pain/opioid response
11. **Rare Diseases** (29) - Rare conditions
12. **Mental Health** (25) - Psychiatric genetics
13. **Dermatology** (37) - Skin and hair
14. **Vision & Hearing** (33) - Sensory genetics
15. **Fertility** (31) - Reproductive health
16. **Nutrition** (34) - Nutrigenomics
17. **Fitness** (30) - Athletic performance
18. **Neurogenetics** (28) - Cognition/behavior
19. **Longevity** (30) - Aging markers
20. **Immunity** (43) - HLA and immune
21. **Ancestry AIMs** (24) - Admixture markers

## Agent Integration

The `agent_summary.json` provides:

```json
{
  "critical_alerts": [],
  "high_priority": [],
  "medium_priority": [],
  "pharmacogenomics_alerts": [],
  "apoe_status": {},
  "polygenic_risk_scores": {},
  "haplogroups": {
    "mtDNA": {"haplogroup": "H", "lineage": "maternal"},
    "Y_DNA": {"haplogroup": "R1b", "lineage": "paternal"}
  },
  "ancestry": {
    "composition": {},
    "admixture": {}
  },
  "hereditary_cancer": {},
  "autoimmune_risk": {},
  "pain_sensitivity": {},
  "lifestyle_recommendations": {
    "diet": [],
    "exercise": [],
    "supplements": [],
    "avoid": []
  },
  "drug_interaction_matrix": {},
  "data_quality": {}
}
```

## Critical Findings (Always Alert User)

### Pharmacogenomics
- **DPYD** variants - 5-FU/capecitabine FATAL toxicity risk
- **HLA-B*5701** - Abacavir hypersensitivity
- **HLA-B*1502** - Carbamazepine SJS (certain populations)
- **MT-RNR1** - Aminoglycoside-induced deafness

### Hereditary Cancer
- **BRCA1/BRCA2** pathogenic - Breast/ovarian cancer syndrome
- **Lynch syndrome** genes - Colorectal/endometrial cancer
- **TP53** pathogenic - Li-Fraumeni syndrome (multi-cancer)

### Disease Risk
- **APOE ε4/ε4** - ~12x Alzheimer's risk
- **Factor V Leiden** - Thrombosis risk, contraceptive implications
- **HLA-B27** - Ankylosing spondylitis susceptibility (OR ~70)

### Carrier Status
- **CFTR** - Cystic fibrosis (1 in 25 Europeans)
- **HBB** - Sickle cell (1 in 12 African Americans)
- **HEXA** - Tay-Sachs (1 in 30 Ashkenazi Jews)

## Usage Examples

### Basic Analysis
```python
from comprehensive_analysis import main
main()  # Uses command line args
```

### Haplogroup Analysis
```python
from markers.haplogroups import analyze_haplogroups
result = analyze_haplogroups(genotypes)
print(result["mtDNA"]["haplogroup"])  # e.g., "H"
```

### Ancestry
```python
from markers.ancestry_composition import get_ancestry_summary
ancestry = get_ancestry_summary(genotypes)
```

### Cancer Panel
```python
from markers.cancer_panel import analyze_cancer_panel
cancer = analyze_cancer_panel(genotypes)
if cancer["pathogenic_variants"]:
    print("ALERT: Pathogenic variants detected")
```

### Generate PDF
```python
from pdf_report import generate_pdf_report
pdf_path = generate_pdf_report(analysis_results)
```

### Export for Genetic Counselor
```python
from exports import generate_genetic_counselor_export
clinical = generate_genetic_counselor_export(results, "clinical.json")
```

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
- Negative hereditary cancer result does NOT rule out cancer syndrome
- Haplogroup resolution limited without WGS

## When to Recommend Genetic Counseling

- Any pathogenic hereditary cancer variant
- APOE ε4/ε4 genotype
- Multiple critical pharmacogenomic findings
- Carrier status with reproduction implications
- High-risk autoimmune HLA types with symptoms
- Results causing significant user distress
