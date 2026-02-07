# Changelog

All notable changes to the Personal Genomics skill will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [4.4.0] - 2026-02-07

### Added

#### 1000 Genomes Population Comparison
- New `references/1000genomes_frequencies.json` - Population frequency data for 20+ markers from 1000 Genomes Phase 3
- New `markers/population_comparison.py` module with:
  - `compare_to_populations()` - Compare genotypes against all 20 reference populations
  - `find_most_similar_populations()` - Rank populations by genetic similarity
  - `get_marker_population_context()` - Get detailed population context per marker
  - `generate_population_comparison_report()` - Text report generation
  - `get_population_comparison_json()` - Structured JSON for dashboard
- 20 reference populations across 5 superpopulations (EUR, AFR, EAS, SAS, AMR)
- Visual bar charts showing genotype frequencies by population
- Transparent methodology - users see actual data, not black-box percentages

#### Ancient DNA Signals
- New `references/ancient_dna_markers.json` - Markers with ancient DNA associations
- New `markers/ancient_ancestry.py` module with:
  - `detect_ancient_signals()` - Find matches to ancient populations
  - `generate_ancient_dna_report()` - Narrative report with historical context
  - `get_ancient_dna_json()` - Structured JSON for dashboard
  - `get_neanderthal_report()` - Detailed Neanderthal introgression analysis
- Coverage of 6 ancient populations:
  - Western Hunter-Gatherers (WHG) - Mesolithic Europe
  - Eastern Hunter-Gatherers (EHG) - Mesolithic Russia
  - Anatolian Neolithic Farmers (ANF) - Early farmers
  - Yamnaya/Steppe Pastoralists - Bronze Age herders
  - Neanderthal introgression markers
  - Denisovan introgression markers
- Timeline visualization showing human migration history
- PMID citations from published ancient DNA studies

#### Enhanced Dashboard
- New "🌍 Population Comparison" section with:
  - Most Similar Populations ranking with similarity scores
  - Per-marker frequency breakdown with visual bars
  - Continental summary view
  - "Why No Percentages?" educational explainer
  - Tab-based navigation between views
- New "🏛️ Ancient DNA" section integrated with existing ancestry
- Updated navigation with new sections
- Added CSS for population bars and ancient DNA timeline

#### Documentation
- New `references/METHODOLOGY.md` explaining:
  - Why we show comparisons instead of percentages
  - 1000 Genomes Project reference
  - Ancient DNA study sources
  - Neanderthal introgression context
  - How to interpret results
  - Limitations of consumer DNA testing

### Changed
- Updated version to 4.4.0 across all files
- Dashboard now shows population comparison alongside ancestry signals
- Agent summary now includes `population_comparison`, `ancient_dna`, and `neanderthal` sections
- Updated SKILL.md with v4.4.0 feature highlights

### Fixed
- Test compatibility for ancestry module API changes
- Backwards compatibility aliases for `ANCESTRY_INFORMATIVE_MARKERS`

## [4.2.0] - 2026-02-07

### Added

#### Interactive Web Dashboard
- New `dashboard/index.html` - self-contained HTML/CSS/JS visualization
- Auto-generated `dashboard.html` with every analysis run
- Sections: Overview, Pharmacogenomics, Health Risks (PRS), Traits, Ancestry, Carrier Status, Sleep Profile, Athletic Profile, UV/Skin, Dietary Interactions
- Dark mode toggle with persistence
- Mobile responsive design
- Collapsible sections for easy navigation
- Search/filter functionality for tables and traits
- Export to PDF via print functionality
- Drag & drop JSON loading for standalone use
- Zero external dependencies - works completely offline
- Inline CSS/JS for single-file distribution

#### Code Quality Improvements
- Complete type hints for all public functions
- TypedDict definitions for complex return types (GenotypeData, MarkerInfo, AnalysisResult, APOEResult, PRSResult)
- Google-style docstrings for all public functions
- Added `py.typed` PEP 561 marker file
- Module-level docstrings explaining purpose

#### Defensive Programming
- Input validation for all public functions
- `validate_rsid()` - validates rsID format
- `validate_genotype()` - validates genotype format
- `validate_filepath()` - validates and resolves file paths
- `sanitize_genotype()` - cleans and normalizes genotype strings
- Graceful handling of missing/malformed genotype data
- Informative error messages with suggested fixes
- Logging for warnings and recoverable issues

#### API Improvements
- New `analyze_dna_file()` function as main entry point for programmatic usage
- New `generate_dashboard()` function for standalone dashboard generation
- `--open` CLI flag to auto-open dashboard in browser
- `--no-dashboard` CLI flag to skip dashboard generation
- Better separation of CLI and library usage

#### Testing
- 30 new edge case tests in `tests/test_edge_cases.py`
- Input validation tests
- Empty input handling tests
- Malformed data tests
- Encoding handling tests (UTF-8, Latin-1)
- APOE edge case tests
- PRS edge case tests
- Report generation tests
- Dashboard generation tests
- Integration tests
- Total: 200 tests, all passing

### Changed
- Updated version to 4.2.0
- Improved error handling throughout codebase
- Better type safety with explicit type annotations
- More defensive code patterns
- Cleaner separation of concerns

### Fixed
- Better handling of edge cases in file loading
- Improved encoding error handling
- Fixed potential issues with missing marker data

## [4.1.0] - 2026-02-06

### Added
- Medication Interaction Checker with brand/generic name support
- Sleep Optimization Profile (chronotype, caffeine metabolism)
- Dietary Interaction Matrix (caffeine, alcohol, lactose, gluten, etc.)
- Athletic Performance Profiling (power vs endurance, recovery, injury risk)
- UV Sensitivity Calculator (skin type, SPF recommendations)
- Natural Language Explanations for all findings
- Telomere Length Estimation
- Research Variant Flagging
- Extended PRS markers
- Extended carrier status markers
- Extended pharmacogenomics markers

## [4.0.0] - 2026-02-06

### Added
- Haplogroup Analysis (mtDNA and Y-DNA)
- Ancestry Composition with population comparisons
- Expanded Hereditary Cancer Panel (BRCA1/2, Lynch syndrome, etc.)
- Autoimmune HLA Associations
- Pain Sensitivity markers
- PDF Report Generation
- Data Quality Metrics
- Integration & Export formats

## [3.0.0] - Initial Release

### Added
- Core genetic analysis with 1600+ markers
- Pharmacogenomics (CPIC Level 1A)
- Polygenic Risk Scores for 10+ conditions
- Carrier Status screening
- Health Risks assessment
- Traits analysis
- Nutrition markers
- Fitness markers
- Neurogenetics markers
- Longevity markers
- Immunity markers
- Rare Diseases panel
- Mental Health markers
- Dermatology markers
- Vision & Hearing markers
- Fertility markers
- Agent-friendly JSON output
- Human-readable text reports
- Multi-format support (23andMe, AncestryDNA, VCF, etc.)
