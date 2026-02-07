# Changelog

All notable changes to the Personal Genomics skill will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [4.4.1] - 2026-02-07

### Fixed

#### Ancient DNA Matching - MAJOR SNP Expansion
- **BREAKING FIX**: Expanded from ~11 total SNPs to 41 unique SNPs across 32 ancient individuals
- Each ancient individual now has 27-41 SNPs (up from 4-9)
- 38 SNPs now overlap with AncestryDNA data (up from 4-5)
- Similarity scores now vary meaningfully (62%-78% range vs. meaningless "100%")

#### Confidence Levels Added
- **High confidence** (🟢): 30+ shared SNPs - strong statistical power
- **Medium confidence** (🟡): 20-29 shared SNPs - good for patterns  
- **Low confidence** (🟠): 10-19 shared SNPs - interpret with caution
- Confidence displayed prominently on all match cards and culture affinities

#### Dashboard Updates
- Match cards now show "X/Y identical" instead of just percentage
- Confidence emoji (🟢🟡🟠🔴) displayed on each card
- Culture affinities show SNP range and sample size (n=X)
- Removed misleading "similarity" label, replaced with actual counts

#### Expanded SNP Categories
- **Pigmentation**: rs1426654, rs16891982, rs12913832, rs1042602, rs1800407, rs12896399, rs7495174, rs4778138, rs1393350, rs6119471, rs12203592, rs2733832, rs1800414
- **Hair/MC1R**: rs1805007, rs1805008, rs1805005, rs1805006, rs2228479, rs11547464, rs885479
- **Diet/metabolism**: rs4988235, rs182549, rs174546
- **Morphology**: rs3827760, rs17822931, rs260690, rs10843090
- **Ancestry-informative**: rs2814778, rs10962599, rs7657799, rs4833103, rs2891
- **Misc phenotype**: rs671, rs1229984, rs4680, rs1799971, rs6265, rs4570625, rs53576

### Changed
- `ancient_individuals.json` completely rebuilt with accurate population-typical genotypes
- `ancient_matching.py` rewritten with confidence scoring system
- Weighted scoring now applies to 11 ancestry-informative SNPs
- Results now ranked by percentile across all 32 ancient individuals

## [4.4.0] - 2026-02-07

### Added

#### Ancient Individual Matching (YourTrueAncestry Clone)
- New `references/ancient_individuals.json` - 30+ ancient genomes from published studies
- New `references/ancient_cultures.json` - 15+ archaeological cultures for affinity matching
- New `markers/ancient_matching.py` module with:
  - `calculate_genetic_distance()` - IBS-based similarity calculation
  - `find_closest_ancients()` - Rank ancient individuals by similarity
  - `match_to_cultures()` - Calculate cultural affinities
  - `generate_ancient_matches_report()` - Detailed text report
  - `get_ancient_matches_json()` - Structured JSON for dashboard
- **Coverage includes:**
  - Mesolithic: Cheddar Man, Loschbour, La Braña
  - Neolithic: Barcın Farmer, Stuttgart LBK, British Neolithic
  - Bronze Age: Ötzi, Yamnaya, Bell Beaker, Corded Ware
  - Iron Age: Hallstatt Celtic, British Iron Age
  - Historical: Roman Britain, Anglo-Saxons, Vikings, Medieval
- **Better than paid services:**
  - ✅ Full transparency - methodology completely open
  - ✅ PMIDs for every sample - verify the science
  - ✅ Downloadable raw data - export anytime
  - ✅ Free forever - no paywall
  - ✅ Open source - fully auditable code

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
