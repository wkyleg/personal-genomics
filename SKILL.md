# Personal Genomics

An OpenClaw skill for local genetic data analysis.

## Description

This skill enables agents to analyze raw DNA files from consumer genetic testing services (23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA). The agent can help users understand health markers, pharmacogenomic variants, ancestry composition, and ancient DNA connections.

All processing runs locally with no network requests. Genetic data never leaves the user's machine.

## Requirements

Python 3.8 or later with pandas and numpy. Optional dependencies include PLINK for population genetics and bcftools for VCF processing.

## Primary Use Cases

**Health inference**: Analyze 800+ markers for disease risk, APOE status, MTHFR variants, and other health-relevant polymorphisms. Results include evidence grades (strong/moderate/preliminary) to help agents calibrate confidence.

**Pharmacogenomics**: Identify drug metabolism variants (CYP2C19, CYP2C9, CYP2D6, VKORC1) affecting response to warfarin, clopidogrel, statins, and other medications.

**Ancestry**: Determine Y-DNA and mtDNA haplogroups across all major lineages. Analyze ancestry-informative markers for European, African, East Asian, South Asian, Middle Eastern, Native American, and Oceanian populations.

**Ancient DNA**: Detect Neanderthal and Denisovan introgression markers. Compare against Allen Ancient DNA Resource when available.

## Output

The skill produces human-readable markdown reports and structured JSON with actionable fields:

```json
{
  "high_priority_alerts": [...],
  "medication_alerts": [...],
  "supplement_considerations": [...],
  "lifestyle_recommendations": [...]
}
```

Priority levels and action types help agents determine what to surface to users.

## Quick Start

```bash
python comprehensive_analysis.py ~/Downloads/23andMe.txt
```

Results written to `~/dna-analysis/reports/`.

## Scripts

| Script | Purpose |
|--------|---------|
| comprehensive_analysis.py | Full analysis with agent-friendly output |
| analyze_dna.py | Core health and pharmacogenomic markers |
| ethnicity_analysis.py | Multi-ethnic ancestry analysis |
| neanderthal_analysis.py | Archaic introgression detection |
| ancient_dna.py | Ancient population markers |
| supplement_protocol.py | Genetically-informed supplement suggestions |
| parental_inference.py | Haplogroup determination |
| convert_to_plink.py | Export to PLINK format |

## Marker Sources

ClinVar (clinical variants), PharmGKB (pharmacogenomics), GWAS Catalog (association studies), and peer-reviewed literature. Citations included with each marker.

## Privacy

Zero network requests. All analysis local. No data collection.

## Limitations

SNP arrays miss rare variants and structural mutations. Results are informational, not diagnostic. Genetic risk is probabilistic. Consult healthcare providers for medical decisions.
