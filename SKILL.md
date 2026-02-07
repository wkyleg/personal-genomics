# Personal Genomics Skill

Comprehensive local DNA analysis with 2000+ markers. Privacy-first genetic analysis for AI agents.

## Quick Start

```bash
python comprehensive_analysis.py /path/to/dna_file.txt
```

## Supported Files

- 23andMe, AncestryDNA, MyHeritage, FTDNA
- VCF files (whole genome/exome)
- Any tab-delimited rsid format

## Output Location

`~/dna-analysis/reports/`

- `agent_summary.json` - AI-optimized, priority-sorted
- `full_analysis.json` - Complete data
- `report.txt` - Human-readable

## Marker Categories

1. **Pharmacogenomics** - Drug metabolism (CYP450, DPYD, TPMT, HLA)
2. **Polygenic Risk Scores** - CAD, T2D, cancer, Alzheimer's, etc.
3. **Carrier Status** - CF, sickle cell, Tay-Sachs, 35+ diseases
4. **Health Risks** - APOE, Factor V Leiden, AMD, celiac
5. **Traits** - Eye color, taste, alcohol flush, chronotype
6. **Nutrition** - MTHFR, vitamin D, caffeine, omega-3
7. **Fitness** - ACTN3, injury risk, recovery
8. **Neurogenetics** - COMT, BDNF, dopamine, serotonin
9. **Longevity** - FOXO3, TERT, aging markers
10. **Immunity** - HLA, autoimmunity, infection risk

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
  "notable_traits": []
}
```

## Key Markers

### Pharmacogenomics (CPIC Level 1A)
- CYP2D6 - Codeine, tramadol, tamoxifen
- CYP2C19 - Clopidogrel (Plavix)
- CYP2C9/VKORC1 - Warfarin dosing
- DPYD - 5-FU (CRITICAL - fatal toxicity risk)
- TPMT/NUDT15 - Azathioprine
- HLA-B*5701 - Abacavir
- HLA-B*1502 - Carbamazepine (Asian ancestry)

### Disease Risk
- APOE - Alzheimer's, cardiovascular
- Factor V Leiden - Thrombosis
- BRCA1/2 proxies - Breast cancer
- 8q24 locus - Multiple cancers

### Carrier Status
- CFTR - Cystic fibrosis (1 in 25 Europeans)
- HBB - Sickle cell (1 in 12 African Americans)
- HEXA - Tay-Sachs (1 in 30 Ashkenazi Jews)

## Privacy

- All analysis runs locally
- Zero network requests
- No data leaves the machine

## Limitations

- Consumer arrays miss rare variants
- Results are probabilistic, not deterministic
- Not a medical diagnosis
- Consult healthcare providers for medical decisions
