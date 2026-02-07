<p align="center">
  <img src="logo.svg" width="200" height="200" alt="Personal Genomics">
</p>

<h1 align="center">Personal Genomics</h1>

<p align="center">
  <strong>Local genetic analysis for OpenClaw agents</strong><br>
  Analyze consumer genetic data with complete privacy
</p>

---

## Overview

Personal Genomics is an OpenClaw skill that enables AI agents to analyze raw genetic data files from consumer testing services. The skill processes genotype data locally, extracting clinically relevant variants, ancestry markers, and pharmacogenomic information that agents can use to help users understand their genetic profile.

The analysis covers approximately 800 well-characterized genetic markers drawn from peer-reviewed literature and curated databases including ClinVar, PharmGKB, and the NHGRI-EBI GWAS Catalog. Each marker includes citations to primary sources, enabling verification of the underlying evidence.

This tool is designed for educational and research purposes. It is not a medical device and should not be used for clinical decision-making without consultation with qualified healthcare professionals.

## Privacy and Data Handling

**All genetic analysis occurs locally on the user's machine.** The skill makes no network requests and transmits no data to external servers. Raw genetic files remain under the user's exclusive control throughout the analysis process.

This design reflects the sensitive nature of genetic information. Unlike many web-based genetic analysis tools, Personal Genomics requires no account creation, collects no usage data, and maintains no connection to cloud services. Users who wish to verify this can audit the source code directly.

Genetic data files should be stored securely and shared cautiously. Even de-identified genetic data can potentially be re-identified through comparison with public databases or relatives' genetic profiles.

## Supported Data Sources

The skill accepts raw data exports from major direct-to-consumer genetic testing providers:

- 23andMe (all chip versions)
- AncestryDNA
- MyHeritage DNA
- FamilyTreeDNA
- LivingDNA
- Nebula Genomics

Most providers offer raw data downloads through account settings. The resulting files contain rsID identifiers, chromosomal positions, and genotype calls in tab-delimited format. The parser handles format variations across providers automatically.

## Technical Implementation

The analysis pipeline builds on established open-source genomics tools:

**PLINK** (Purcell et al., 2007) provides population genetics calculations, linkage disequilibrium analysis, and format conversions. The skill includes utilities for converting consumer genetic data into PLINK binary format, enabling integration with the broader genomics software ecosystem.

**bcftools** from the Samtools project handles VCF file operations when working with sequencing data or annotated variant files.

Core data processing uses **pandas** and **numpy**, which efficiently handle the 600,000 to 700,000 SNPs typically present in consumer genotyping arrays.

## Analysis Modules

### Health-Associated Variants

The skill evaluates markers with established associations to disease risk, including variants implicated in cardiovascular disease, metabolic disorders, cancer predisposition, and autoimmune conditions. APOE genotype, relevant to Alzheimer's disease and lipid metabolism, is inferred from the rs429358 and rs7412 polymorphisms.

Results include effect sizes and evidence classifications (strong, moderate, or preliminary) based on the quality and replication status of underlying studies. Agents can use these classifications to calibrate how findings are presented to users.

### Pharmacogenomics

Genetic variation substantially influences drug metabolism and response. The skill identifies variants affecting cytochrome P450 enzymes (CYP2C19, CYP2C9, CYP2D6, CYP1A2), drug transporters (SLCO1B1), and drug targets (VKORC1).

Clinically actionable variants include those affecting warfarin dosing, clopidogrel efficacy, statin-induced myopathy risk, and SSRI metabolism. The Clinical Pharmacogenetics Implementation Consortium (CPIC) guidelines inform interpretation of these variants.

Users should share pharmacogenomic findings with prescribers, particularly before initiating medications with known genetic interactions.

### Ancestry and Population Genetics

The skill determines Y-chromosome and mitochondrial DNA haplogroups, which trace paternal and maternal lineages respectively. Coverage spans all major haplogroup clades found in global populations, from African-origin lineages (Y-DNA A, B; mtDNA L0-L6) through more recent branches specific to European, Asian, Oceanian, and American populations.

Ancestry-informative markers enable estimation of continental ancestry proportions. The reference populations include European, Sub-Saharan African, East Asian, South Asian, Middle Eastern, Native American, and Oceanian groups.

### Archaic Introgression

Modern humans outside Africa carry approximately 1-2% Neanderthal DNA from interbreeding events roughly 50,000-60,000 years ago. The skill identifies known introgression markers, including variants affecting immune function, skin and hair characteristics, and metabolism.

For users interested in deep ancestry, the skill can compare results against the Allen Ancient DNA Resource (AADR) database when that reference dataset is available locally.

## Output Format

Analysis produces both human-readable markdown reports and structured JSON designed for agent consumption. The JSON format includes actionable metadata:

```json
{
  "analysis_date": "2024-01-15T14:30:00",
  "total_markers_analyzed": 847,
  "high_priority_alerts": [
    {
      "gene": "CYP2C19",
      "rsid": "rs4244285",
      "genotype": "AG",
      "evidence_level": "strong",
      "action_type": "medical_alert",
      "recommendations": [
        "Discuss with prescriber before taking clopidogrel",
        "Alternative antiplatelet therapy may be indicated"
      ],
      "references": ["PMID:21716271", "CPIC guideline"]
    }
  ],
  "medication_alerts": [...],
  "supplement_considerations": [...],
  "lifestyle_recommendations": [...]
}
```

Priority levels (high, medium, low, informational) and action types (medical_alert, screening, supplementation, lifestyle_modification) help agents determine what information to emphasize.

## Usage

Run comprehensive analysis:

```
python comprehensive_analysis.py /path/to/genetic_data.txt
```

Individual analysis modules:

```
python analyze_dna.py <file>              # Health markers
python pharmacogenomics.py <file>         # Drug metabolism
python ethnicity_analysis.py <file>       # Ancestry composition
python neanderthal_analysis.py <file>     # Archaic introgression
python supplement_protocol.py <file>      # Nutritional genetics
```

Output files are written to `~/dna-analysis/reports/`.

## Limitations and Disclaimers

**This software is not a medical device.** Results are intended for educational and research purposes only. Genetic risk predictions are probabilistic, not deterministic. A risk-associated variant indicates elevated probability relative to population baseline, not certainty of disease.

Consumer genotyping arrays have inherent limitations. They assess only a subset of genetic variation (typically 0.02% of the genome) and miss rare variants, structural variants, and many clinically significant mutations that would be detected by clinical-grade sequencing.

Genetic associations identified in research populations may not generalize across all ancestries. Many variants in current databases were characterized primarily in European-descent populations and may have different effects or frequencies in other groups.

Environmental factors, lifestyle, epigenetics, and gene-gene interactions all influence phenotypic outcomes. Genetic information should be interpreted as one component of a comprehensive health assessment, not as a standalone predictor.

Users considering medical decisions based on genetic information should consult with genetic counselors or healthcare providers qualified to interpret results in clinical context.

## Contributing

Contributions are welcome. When adding markers, please include:

- Primary literature citations (PMID preferred)
- ClinVar accession numbers where available
- Effect allele, risk direction, and effect size
- Population(s) in which the association was characterized
- Any known ancestry-specific effects

Test contributions with genetic data from diverse ancestry backgrounds to ensure broad applicability.

## License

MIT License. See LICENSE file.

The software includes additional terms regarding genetic analysis: results are provided for informational purposes only and should not be used as the basis for medical decisions without professional consultation.

## References

Purcell S, et al. (2007). PLINK: A tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics, 81(3), 559-575.

Landrum MJ, et al. (2018). ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Research, 46(D1), D1062-D1067.

Whirl-Carrillo M, et al. (2021). An evidence-based framework for evaluating pharmacogenomics knowledge for personalized medicine. Clinical Pharmacology & Therapeutics, 110(3), 563-572.

Buniello A, et al. (2019). The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics. Nucleic Acids Research, 47(D1), D1005-D1012.
