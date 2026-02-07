"""
Haplogroup Analysis Module
Mitochondrial DNA (mtDNA) and Y-Chromosome Haplogroups

Sources:
- PhyloTree (mtDNA tree) - https://www.phylotree.org/
- ISOGG Y-DNA Haplogroup Tree - https://isogg.org/tree/
- Published phylogenetic studies

╔══════════════════════════════════════════════════════════════════════════════╗
║  ⚠️  IMPORTANT LIMITATION - CONSUMER ARRAY DATA                              ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Consumer DNA arrays (23andMe, AncestryDNA, etc.) capture only a SMALL       ║
║  SUBSET of haplogroup-defining markers. These results are:                   ║
║                                                                              ║
║  • LOW CONFIDENCE - indicative only, not definitive                          ║
║  • LIMITED RESOLUTION - cannot determine deep subclades                       ║
║  • POTENTIALLY INCORRECT - missing markers can lead to misclassification     ║
║                                                                              ║
║  For accurate haplogroup determination, use dedicated testing:               ║
║  • Y-DNA: FTDNA Y-37/Y-111/Big Y, YFull, Nebula WGS                         ║
║  • mtDNA: FTDNA mtDNA Full Sequence, YFull, Nebula WGS                      ║
╚══════════════════════════════════════════════════════════════════════════════╝

PMIDs:
- van Oven M, Kayser M. 2009. Updated comprehensive phylogenetic tree of global 
  human mtDNA variation. PMID: 19165223
- Karafet TM et al. 2008. New binary polymorphisms reshape Y-chromosome haplogroup 
  tree. PMID: 18285812
- Poznik GD et al. 2016. Punctuated bursts in human male demography. PMID: 27654910
"""

from typing import Dict, List, Optional, Any, Tuple

# Disclaimer text for all haplogroup output
HAPLOGROUP_DISCLAIMER = """
⚠️ CONSUMER ARRAY LIMITATION - LOW CONFIDENCE RESULT

This haplogroup call is based on LIMITED markers from a consumer DNA array.
Consumer arrays cannot reliably determine haplogroups because:

1. They test only ~20-50 haplogroup markers (vs thousands defining the tree)
2. Missing branch-defining SNPs can cause misclassification
3. Resolution is limited to major haplogroups only (no subclades)

CONFIDENCE LEVEL: LOW - INDICATIVE ONLY

For accurate, verified haplogroup determination:
• Y-DNA: FTDNA Big Y-700, YFull WGS analysis
• mtDNA: FTDNA mtDNA Full Sequence, YFull mtDNA
• Both: Whole genome sequencing (Nebula, Dante Labs)

Do not rely on these results for genealogical research or family matching.
"""

# =============================================================================
# MITOCHONDRIAL DNA (mtDNA) HAPLOGROUPS - MATERNAL LINEAGE
# =============================================================================

MTDNA_MARKERS = {
    # Major defining SNPs for mtDNA haplogroups
    # Format: rsid -> {haplogroup, position, ref, alt, branch_defining}
    # PMIDs: 19165223 (PhyloTree), 18498415 (mtDNA phylogeny)
    
    # Haplogroup L (African origin - all humans descend from)
    "rs2853499": {
        "position": "m.263A>G",
        "haplogroup": "L3",
        "branch": "L3 and descendants",
        "alleles": {"G": "L3+", "A": "root"},
        "significance": "ancestral_branch",
        "pmid": ["19165223"]
    },
    
    # Haplogroup M - Out of Africa (Asia, Americas)
    "rs28358571": {
        "position": "m.10398A>G",
        "haplogroup": "M/N indicator",
        "branch": "Macrohaplogroup split",
        "alleles": {"G": "M_branch", "A": "N_branch"},
        "note": "Major split in human migration",
        "pmid": ["19165223", "18498415"]
    },
    
    # Haplogroup N - Out of Africa (Europe, West Asia, partial Americas)
    "rs2857284": {
        "position": "m.8701A>G",
        "haplogroup": "N",
        "branch": "Macrohaplogroup N",
        "alleles": {"G": "N+"},
        "pmid": ["19165223"]
    },
    
    # Haplogroup R - Derived from N (Europe, South Asia)
    "rs2853495": {
        "position": "m.12705C>T",
        "haplogroup": "R",
        "branch": "Haplogroup R",
        "alleles": {"T": "R+"},
        "pmid": ["19165223"]
    },
    
    # Haplogroup H - Most common European (~40-50%)
    "rs2853511": {
        "position": "m.7028C>T",
        "haplogroup": "H",
        "branch": "Haplogroup H",
        "alleles": {"C": "H+", "T": "non-H"},
        "frequency": {"Europe": 0.45, "Middle East": 0.20},
        "pmid": ["19165223", "12730389"]
    },
    "rs2853512": {
        "position": "m.2706A>G",
        "haplogroup": "H",
        "branch": "H confirmation",
        "alleles": {"A": "H+"},
        "pmid": ["19165223"]
    },
    
    # Haplogroup V - European (Basque, Scandinavian)
    "rs2857291": {
        "position": "m.4580G>A",
        "haplogroup": "V",
        "branch": "Haplogroup V",
        "alleles": {"A": "V+"},
        "frequency": {"Basque": 0.12, "Scandinavia": 0.06},
        "pmid": ["19165223", "10839820"]
    },
    
    # Haplogroup U - European, West Asian
    "rs41323649": {
        "position": "m.11467A>G",
        "haplogroup": "U",
        "branch": "Haplogroup U",
        "alleles": {"G": "U+"},
        "pmid": ["19165223"]
    },
    "rs28359178": {
        "position": "m.12308A>G",
        "haplogroup": "U",
        "branch": "U confirmation",
        "alleles": {"G": "U+"},
        "pmid": ["19165223"]
    },
    
    # Haplogroup K - European (Ashkenazi Jewish enriched)
    "rs28357370": {
        "position": "m.9055G>A",
        "haplogroup": "K",
        "branch": "Haplogroup K (U8b)",
        "alleles": {"A": "K+"},
        "frequency": {"Ashkenazi": 0.32, "Europe": 0.06},
        "pmid": ["19165223", "16689639"]
    },
    
    # Haplogroup J - European, Middle Eastern
    "rs2853504": {
        "position": "m.295C>T",
        "haplogroup": "J",
        "branch": "Haplogroup J",
        "alleles": {"T": "J+"},
        "frequency": {"Europe": 0.09, "Middle East": 0.12},
        "pmid": ["19165223"]
    },
    "rs2853505": {
        "position": "m.489T>C",
        "haplogroup": "J",
        "branch": "J confirmation",
        "alleles": {"C": "J+"},
        "pmid": ["19165223"]
    },
    
    # Haplogroup T - European, Middle Eastern
    "rs2853506": {
        "position": "m.709G>A",
        "haplogroup": "T",
        "branch": "Haplogroup T",
        "alleles": {"A": "T+"},
        "frequency": {"Europe": 0.09},
        "pmid": ["19165223"]
    },
    
    # Haplogroup I - European (Scandinavian)
    "rs28358569": {
        "position": "m.10034T>C",
        "haplogroup": "I",
        "branch": "Haplogroup I",
        "alleles": {"C": "I+"},
        "frequency": {"Scandinavia": 0.04},
        "pmid": ["19165223"]
    },
    
    # Haplogroup W - European, West Asian
    "rs2853507": {
        "position": "m.1243T>C",
        "haplogroup": "W",
        "branch": "Haplogroup W",
        "alleles": {"C": "W+"},
        "pmid": ["19165223"]
    },
    
    # Haplogroup X - Native American (minor), European
    "rs2853508": {
        "position": "m.1719G>A",
        "haplogroup": "X",
        "branch": "Haplogroup X",
        "alleles": {"A": "X+"},
        "frequency": {"Native American": 0.03, "Europe": 0.02},
        "pmid": ["19165223", "12547090"]
    },
    
    # Haplogroup A - East Asian, Native American
    "rs2853510": {
        "position": "m.663A>G",
        "haplogroup": "A",
        "branch": "Haplogroup A",
        "alleles": {"G": "A+"},
        "frequency": {"East Asia": 0.10, "Native American": 0.25},
        "pmid": ["19165223"]
    },
    
    # Haplogroup B - East Asian, Polynesian, Native American
    "rs2853509": {
        "position": "m.8281-8289del",
        "haplogroup": "B",
        "branch": "Haplogroup B",
        "alleles": {"del": "B+"},
        "frequency": {"Polynesia": 0.90, "Native American": 0.20},
        "pmid": ["19165223"]
    },
    
    # Haplogroup C - East Asian, Native American
    "rs2853496": {
        "position": "m.14318T>C",
        "haplogroup": "C",
        "branch": "Haplogroup C",
        "alleles": {"C": "C+"},
        "frequency": {"Native American": 0.30, "East Asia": 0.15},
        "pmid": ["19165223"]
    },
    
    # Haplogroup D - East Asian, Native American
    "rs2853497": {
        "position": "m.5178C>A",
        "haplogroup": "D",
        "branch": "Haplogroup D",
        "alleles": {"A": "D+"},
        "frequency": {"East Asia": 0.25, "Native American": 0.20},
        "pmid": ["19165223"]
    },
    
    # Haplogroup F - East Asian
    "rs2853498": {
        "position": "m.10310G>A",
        "haplogroup": "F",
        "branch": "Haplogroup F",
        "alleles": {"A": "F+"},
        "frequency": {"Southeast Asia": 0.20},
        "pmid": ["19165223"]
    },
    
    # African haplogroups
    # L0 - Oldest human lineage (Khoisan)
    "rs3928306": {
        "position": "m.3594C>T",
        "haplogroup": "L0",
        "branch": "Haplogroup L0",
        "alleles": {"T": "L0+"},
        "frequency": {"Khoisan": 0.60},
        "pmid": ["19165223", "20705718"]
    },
    
    # L1 - African
    "rs3135027": {
        "position": "m.769G>A",
        "haplogroup": "L1",
        "branch": "Haplogroup L1",
        "alleles": {"A": "L1+"},
        "frequency": {"Central Africa": 0.30},
        "pmid": ["19165223"]
    },
    
    # L2 - Most common African haplogroup
    "rs3135028": {
        "position": "m.7175T>C",
        "haplogroup": "L2",
        "branch": "Haplogroup L2",
        "alleles": {"C": "L2+"},
        "frequency": {"West Africa": 0.35, "African American": 0.30},
        "pmid": ["19165223"]
    },
}

# =============================================================================
# Y-CHROMOSOME HAPLOGROUPS - PATERNAL LINEAGE (Males only)
# =============================================================================

YCHROMOSOME_MARKERS = {
    # Major Y-DNA haplogroup defining SNPs
    # These require Y-chromosome coverage (males only)
    # PMIDs: 18285812 (Y tree), 27654910 (Y phylogeny)
    
    # Haplogroup A - Oldest Y lineage (African)
    "rs2032658": {
        "gene": "Y",
        "position": "M91",
        "haplogroup": "A",
        "branch": "Root haplogroup A",
        "alleles": {"T": "A+"},
        "sex": "male",
        "origin": "Africa",
        "age_kya": 270,
        "pmid": ["18285812"]
    },
    
    # Haplogroup B - African
    "rs9786184": {
        "gene": "Y",
        "position": "M60",
        "haplogroup": "B",
        "branch": "Haplogroup B",
        "alleles": {"T": "B+"},
        "sex": "male",
        "origin": "Africa",
        "age_kya": 75,
        "pmid": ["18285812"]
    },
    
    # Haplogroup E - African, spread to Mediterranean
    "rs2032604": {
        "gene": "Y",
        "position": "M96",
        "haplogroup": "E",
        "branch": "Haplogroup E",
        "alleles": {"G": "E+"},
        "sex": "male",
        "origin": "Africa",
        "age_kya": 65,
        "frequency": {"Africa": 0.60, "Middle East": 0.20, "Europe": 0.10},
        "pmid": ["18285812", "21481273"]
    },
    "rs9306841": {
        "gene": "Y",
        "position": "M35",
        "haplogroup": "E1b1b",
        "branch": "E1b1b (E-M35)",
        "alleles": {"C": "E1b1b+"},
        "sex": "male",
        "note": "Common in Mediterranean, Near East",
        "pmid": ["18285812"]
    },
    
    # Haplogroup C - Asia, Oceania, Americas
    "rs35284970": {
        "gene": "Y",
        "position": "M130",
        "haplogroup": "C",
        "branch": "Haplogroup C",
        "alleles": {"T": "C+"},
        "sex": "male",
        "origin": "Asia",
        "age_kya": 60,
        "frequency": {"Mongolia": 0.50, "Australia": 0.60},
        "pmid": ["18285812"]
    },
    
    # Haplogroup D - East Asia (Tibetan, Japanese)
    "rs2032597": {
        "gene": "Y",
        "position": "M174",
        "haplogroup": "D",
        "branch": "Haplogroup D",
        "alleles": {"T": "D+"},
        "sex": "male",
        "origin": "Asia",
        "frequency": {"Tibet": 0.50, "Japan": 0.35},
        "pmid": ["18285812"]
    },
    
    # Haplogroup G - Caucasus, Mediterranean
    "rs2032636": {
        "gene": "Y",
        "position": "M201",
        "haplogroup": "G",
        "branch": "Haplogroup G",
        "alleles": {"T": "G+"},
        "sex": "male",
        "origin": "Caucasus",
        "age_kya": 50,
        "frequency": {"Georgia": 0.60, "Italy": 0.10},
        "pmid": ["18285812", "22232597"]
    },
    
    # Haplogroup I - European (especially Scandinavian)
    "rs2032637": {
        "gene": "Y",
        "position": "M170",
        "haplogroup": "I",
        "branch": "Haplogroup I",
        "alleles": {"A": "I+"},
        "sex": "male",
        "origin": "Europe",
        "age_kya": 40,
        "frequency": {"Scandinavia": 0.40, "Balkans": 0.40},
        "pmid": ["18285812", "21058503"]
    },
    "rs17316180": {
        "gene": "Y",
        "position": "M253",
        "haplogroup": "I1",
        "branch": "I1 (Scandinavian)",
        "alleles": {"C": "I1+"},
        "sex": "male",
        "frequency": {"Scandinavia": 0.35, "Germany": 0.20},
        "pmid": ["18285812"]
    },
    "rs9341302": {
        "gene": "Y",
        "position": "M423",
        "haplogroup": "I2",
        "branch": "I2 (Balkan/Sardinian)",
        "alleles": {"A": "I2+"},
        "sex": "male",
        "frequency": {"Balkans": 0.40, "Sardinia": 0.40},
        "pmid": ["18285812"]
    },
    
    # Haplogroup J - Middle East, Mediterranean, Jewish
    "rs2032638": {
        "gene": "Y",
        "position": "M304",
        "haplogroup": "J",
        "branch": "Haplogroup J",
        "alleles": {"C": "J+"},
        "sex": "male",
        "origin": "Middle East",
        "age_kya": 45,
        "frequency": {"Middle East": 0.40, "Ashkenazi": 0.25},
        "pmid": ["18285812", "25079122"]
    },
    "rs9341296": {
        "gene": "Y",
        "position": "M267",
        "haplogroup": "J1",
        "branch": "J1 (Arabian/Jewish Cohen)",
        "alleles": {"T": "J1+"},
        "sex": "male",
        "frequency": {"Arabia": 0.40, "Jewish": 0.15},
        "pmid": ["18285812"]
    },
    "rs2032639": {
        "gene": "Y",
        "position": "M172",
        "haplogroup": "J2",
        "branch": "J2 (Mediterranean/Fertile Crescent)",
        "alleles": {"T": "J2+"},
        "sex": "male",
        "frequency": {"Greece": 0.25, "Italy": 0.20},
        "pmid": ["18285812"]
    },
    
    # Haplogroup N - Uralic, Siberian, Baltic
    "rs9341278": {
        "gene": "Y",
        "position": "M231",
        "haplogroup": "N",
        "branch": "Haplogroup N",
        "alleles": {"A": "N+"},
        "sex": "male",
        "origin": "Siberia",
        "frequency": {"Finland": 0.60, "Yakutia": 0.90},
        "pmid": ["18285812"]
    },
    
    # Haplogroup O - East Asian
    "rs2032640": {
        "gene": "Y",
        "position": "M175",
        "haplogroup": "O",
        "branch": "Haplogroup O",
        "alleles": {"T": "O+"},
        "sex": "male",
        "origin": "East Asia",
        "age_kya": 35,
        "frequency": {"China": 0.75, "Japan": 0.55},
        "pmid": ["18285812"]
    },
    
    # Haplogroup Q - Native American, Siberian
    "rs2032651": {
        "gene": "Y",
        "position": "M242",
        "haplogroup": "Q",
        "branch": "Haplogroup Q",
        "alleles": {"T": "Q+"},
        "sex": "male",
        "origin": "Central Asia",
        "frequency": {"Native American": 0.90, "Siberia": 0.25},
        "pmid": ["18285812", "12740896"]
    },
    
    # Haplogroup R - Europe, South Asia
    "rs2032641": {
        "gene": "Y",
        "position": "M207",
        "haplogroup": "R",
        "branch": "Haplogroup R",
        "alleles": {"A": "R+"},
        "sex": "male",
        "origin": "Central Asia",
        "age_kya": 30,
        "pmid": ["18285812", "25941402"]
    },
    "rs9786139": {
        "gene": "Y",
        "position": "M17/M198",
        "haplogroup": "R1a",
        "branch": "R1a (Indo-European spread)",
        "alleles": {"C": "R1a+"},
        "sex": "male",
        "frequency": {"Poland": 0.55, "India": 0.35, "Russia": 0.45},
        "pmid": ["18285812", "25941402"]
    },
    "rs9786076": {
        "gene": "Y",
        "position": "M269",
        "haplogroup": "R1b",
        "branch": "R1b (Western European)",
        "alleles": {"C": "R1b+"},
        "sex": "male",
        "frequency": {"Ireland": 0.80, "Spain": 0.65, "France": 0.55},
        "pmid": ["18285812", "25941402"]
    },
    
    # Haplogroup T - Mediterranean, East African
    "rs34126399": {
        "gene": "Y",
        "position": "M184",
        "haplogroup": "T",
        "branch": "Haplogroup T",
        "alleles": {"T": "T+"},
        "sex": "male",
        "frequency": {"Somalia": 0.10, "Fulani": 0.18},
        "pmid": ["18285812"]
    },
}

# =============================================================================
# HAPLOGROUP MIGRATION HISTORY
# =============================================================================

HAPLOGROUP_HISTORY = {
    # mtDNA Haplogroups
    "L0": {
        "type": "mtDNA",
        "origin": "East Africa",
        "age_kya": 150,
        "migration": [
            "Oldest human mitochondrial lineage",
            "Prevalent among Khoisan peoples of southern Africa",
            "Represents deepest branch of human maternal ancestry",
            "Did not participate in Out of Africa migration"
        ],
        "modern_distribution": "Southern Africa (Khoisan), East Africa",
        "pmid": ["19165223", "20705718"]
    },
    "L1": {
        "type": "mtDNA",
        "origin": "Central/West Africa",
        "age_kya": 130,
        "migration": [
            "Ancient African lineage",
            "Common in Central African populations",
            "Some subclades associated with Bantu expansion"
        ],
        "modern_distribution": "Central and West Africa",
        "pmid": ["19165223"]
    },
    "L2": {
        "type": "mtDNA",
        "origin": "West Africa",
        "age_kya": 90,
        "migration": [
            "Most common African mtDNA haplogroup",
            "Expanded with Bantu migrations",
            "High frequency in African diaspora due to slave trade"
        ],
        "modern_distribution": "West Africa, African Americans",
        "pmid": ["19165223"]
    },
    "L3": {
        "type": "mtDNA",
        "origin": "East Africa",
        "age_kya": 70,
        "migration": [
            "Ancestor of all non-African mtDNA",
            "Source of the Out of Africa migration ~60-70 kya",
            "Gave rise to haplogroups M and N"
        ],
        "modern_distribution": "East Africa, and via M/N descendants worldwide",
        "pmid": ["19165223", "17978185"]
    },
    "M": {
        "type": "mtDNA",
        "origin": "East Africa → Asia",
        "age_kya": 60,
        "migration": [
            "Southern route Out of Africa",
            "Rapid coastal migration to South Asia and beyond",
            "Dominant in South Asia, East Asia"
        ],
        "modern_distribution": "South Asia, East Asia, Oceania, Native Americans",
        "pmid": ["19165223"]
    },
    "N": {
        "type": "mtDNA",
        "origin": "East Africa → Middle East",
        "age_kya": 60,
        "migration": [
            "Northern route Out of Africa",
            "Spread through Middle East to Europe and Asia",
            "Gave rise to European haplogroups (H, V, J, T, U, K, etc.)"
        ],
        "modern_distribution": "Europe, Middle East, partial Asia and Americas",
        "pmid": ["19165223"]
    },
    "H": {
        "type": "mtDNA",
        "origin": "Near East → Europe",
        "age_kya": 25,
        "migration": [
            "Re-expanded into Europe after Last Glacial Maximum",
            "Associated with population recovery from Ice Age refugia",
            "Most common European mtDNA haplogroup (40-50%)"
        ],
        "modern_distribution": "Europe (especially Western), Middle East",
        "pmid": ["19165223", "12730389"]
    },
    "U": {
        "type": "mtDNA",
        "origin": "West Asia → Europe",
        "age_kya": 50,
        "migration": [
            "One of earliest European lineages",
            "Present in European hunter-gatherers before farming",
            "Diverse subclades (U5 oldest European-specific)"
        ],
        "modern_distribution": "Europe, South Asia, North Africa",
        "pmid": ["19165223"]
    },
    "K": {
        "type": "mtDNA",
        "origin": "Near East → Europe",
        "age_kya": 30,
        "migration": [
            "Subclade of U8",
            "Enriched in Ashkenazi Jewish populations",
            "Also associated with Neolithic farmers"
        ],
        "modern_distribution": "Europe, especially Ashkenazi Jews",
        "pmid": ["19165223", "16689639"]
    },
    "A": {
        "type": "mtDNA",
        "origin": "East Asia → Americas",
        "age_kya": 25,
        "migration": [
            "One of five founding Native American lineages",
            "Crossed Beringia during Ice Age",
            "Also common in Northeast Asia"
        ],
        "modern_distribution": "Native Americans, East Asia (especially northeast)",
        "pmid": ["19165223"]
    },
    "B": {
        "type": "mtDNA",
        "origin": "East Asia → Pacific/Americas",
        "age_kya": 40,
        "migration": [
            "Dominant Polynesian mtDNA lineage",
            "Part of Austronesian expansion across Pacific",
            "Also founding Native American lineage"
        ],
        "modern_distribution": "Polynesia (very high), Southeast Asia, Native Americans",
        "pmid": ["19165223"]
    },
    
    # Y-DNA Haplogroups
    "R1b": {
        "type": "Y-DNA",
        "origin": "Central Asia → Europe",
        "age_kya": 18,
        "migration": [
            "Associated with Indo-European/Yamnaya expansion",
            "Spread across Europe in Bronze Age",
            "Replaced most prior European Y-DNA lineages",
            "Now dominant in Western Europe"
        ],
        "modern_distribution": "Western Europe (Ireland 80%, Spain 65%)",
        "pmid": ["18285812", "25941402"]
    },
    "R1a": {
        "type": "Y-DNA",
        "origin": "Pontic Steppe → Europe/South Asia",
        "age_kya": 15,
        "migration": [
            "Indo-European expansion, possibly Indo-Iranians",
            "Spread with pastoralist migrations",
            "High in Eastern Europe and South Asia"
        ],
        "modern_distribution": "Eastern Europe, South Asia, Central Asia",
        "pmid": ["18285812", "25941402"]
    },
    "I": {
        "type": "Y-DNA",
        "origin": "Europe (indigenous)",
        "age_kya": 40,
        "migration": [
            "Pre-farming European lineage",
            "Survived in Mesolithic hunter-gatherers",
            "I1 (Scandinavian) and I2 (Balkan) major subclades"
        ],
        "modern_distribution": "Europe (Scandinavia, Balkans, Germany)",
        "pmid": ["18285812", "21058503"]
    },
    "E1b1b": {
        "type": "Y-DNA",
        "origin": "East Africa → Mediterranean",
        "age_kya": 25,
        "migration": [
            "Spread from Africa to Near East and Mediterranean",
            "Associated with Neolithic farmers in some regions",
            "Common in Middle East, North Africa, Mediterranean"
        ],
        "modern_distribution": "Africa, Middle East, Mediterranean Europe",
        "pmid": ["18285812", "21481273"]
    },
    "J": {
        "type": "Y-DNA",
        "origin": "Middle East/Fertile Crescent",
        "age_kya": 45,
        "migration": [
            "Spread with Neolithic farming expansion",
            "J1 associated with Semitic speakers",
            "J2 spread with Mediterranean civilizations"
        ],
        "modern_distribution": "Middle East, Mediterranean, Jewish diaspora",
        "pmid": ["18285812", "25079122"]
    },
    "O": {
        "type": "Y-DNA",
        "origin": "East Asia",
        "age_kya": 35,
        "migration": [
            "Dominant East Asian Y-DNA haplogroup",
            "Expanded with rice agriculture",
            "O2 associated with Han Chinese expansion"
        ],
        "modern_distribution": "East Asia, Southeast Asia",
        "pmid": ["18285812"]
    },
    "Q": {
        "type": "Y-DNA",
        "origin": "Central Asia → Americas",
        "age_kya": 20,
        "migration": [
            "Primary Native American paternal lineage",
            "Crossed Beringia to populate Americas",
            "Also present in Siberia, Central Asia"
        ],
        "modern_distribution": "Native Americans (~90%), Siberia",
        "pmid": ["18285812", "12740896"]
    },
    "N": {
        "type": "Y-DNA",
        "origin": "East Asia → Siberia/Europe",
        "age_kya": 20,
        "migration": [
            "Spread across Siberia and into Northern Europe",
            "Associated with Uralic language family",
            "Dominant in Finland and Baltic peoples"
        ],
        "modern_distribution": "Finland, Baltic, Siberia",
        "pmid": ["18285812"]
    },
}


def determine_mtdna_haplogroup(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Determine mitochondrial DNA haplogroup from available SNPs.
    
    ⚠️ LOW CONFIDENCE - Consumer arrays cannot reliably determine haplogroups.
    Results are INDICATIVE ONLY. Use dedicated mtDNA testing for accuracy.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Dict with haplogroup prediction and confidence
    """
    results = {
        "haplogroup": "Unknown",
        "confidence": "low - consumer array limitation",
        "confidence_note": "Consumer arrays test only ~20 of thousands of haplogroup-defining markers",
        "markers_found": 0,
        "markers_checked": len(MTDNA_MARKERS),
        "branch_evidence": [],
        "possible_subclades": [],
        "history": None,
        "disclaimer": HAPLOGROUP_DISCLAIMER,
        "recommendation": "For accurate mtDNA haplogroup: FTDNA mtDNA Full Sequence or YFull",
        "pmid": ["19165223"]
    }
    
    # Check each marker
    branch_hits = {}
    for rsid, info in MTDNA_MARKERS.items():
        geno = genotypes.get(rsid)
        if not geno:
            continue
            
        results["markers_found"] += 1
        alleles = info.get("alleles", {})
        
        # Check if any allele in genotype matches a haplogroup indicator
        for allele in geno:
            if allele in alleles:
                hg = info.get("haplogroup")
                branch = alleles[allele]
                
                if branch != "root" and "+" in branch:
                    if hg not in branch_hits:
                        branch_hits[hg] = []
                    branch_hits[hg].append({
                        "rsid": rsid,
                        "position": info.get("position"),
                        "genotype": geno,
                        "branch": branch,
                        "pmid": info.get("pmid", [])
                    })
    
    # Determine most likely haplogroup
    if branch_hits:
        # Sort by number of supporting markers
        sorted_hgs = sorted(branch_hits.items(), key=lambda x: len(x[1]), reverse=True)
        top_hg = sorted_hgs[0][0]
        
        results["haplogroup"] = top_hg
        results["branch_evidence"] = sorted_hgs[0][1]
        
        # ALWAYS low confidence for consumer arrays
        # Even with multiple markers, we're missing most defining SNPs
        marker_count = len(sorted_hgs[0][1])
        if marker_count >= 3:
            results["confidence"] = "low - possibly correct major haplogroup"
        elif marker_count == 2:
            results["confidence"] = "very low - limited evidence"
        else:
            results["confidence"] = "very low - single marker only"
        
        # Add possible subclades
        if len(sorted_hgs) > 1:
            results["possible_subclades"] = [hg for hg, _ in sorted_hgs[1:4]]
        
        # Add history
        if top_hg in HAPLOGROUP_HISTORY:
            results["history"] = HAPLOGROUP_HISTORY[top_hg]
    
    return results


def determine_y_haplogroup(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Determine Y-chromosome haplogroup from available SNPs.
    Only applicable for males (requires Y-chromosome data).
    
    ⚠️ LOW CONFIDENCE - Consumer arrays cannot reliably determine haplogroups.
    Results are INDICATIVE ONLY. Use dedicated Y-DNA testing for accuracy.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Dict with haplogroup prediction and confidence
    """
    results = {
        "haplogroup": "Unknown",
        "confidence": "low - consumer array limitation",
        "confidence_note": "Consumer arrays test only ~25 of thousands of Y-DNA SNPs",
        "markers_found": 0,
        "markers_checked": len(YCHROMOSOME_MARKERS),
        "branch_evidence": [],
        "possible_subclades": [],
        "history": None,
        "sex_determination": "unknown",
        "disclaimer": HAPLOGROUP_DISCLAIMER,
        "recommendation": "For accurate Y-DNA haplogroup: FTDNA Big Y-700 or YFull",
        "pmid": ["18285812", "27654910"]
    }
    
    # Check for Y-chromosome presence (indicates male)
    y_markers_found = 0
    branch_hits = {}
    
    for rsid, info in YCHROMOSOME_MARKERS.items():
        geno = genotypes.get(rsid)
        if not geno:
            continue
            
        y_markers_found += 1
        results["markers_found"] += 1
        
        alleles = info.get("alleles", {})
        
        for allele in geno:
            if allele in alleles:
                hg = info.get("haplogroup")
                branch = alleles[allele]
                
                if "+" in branch:
                    if hg not in branch_hits:
                        branch_hits[hg] = []
                    branch_hits[hg].append({
                        "rsid": rsid,
                        "position": info.get("position"),
                        "genotype": geno,
                        "branch": branch,
                        "pmid": info.get("pmid", [])
                    })
    
    # Determine sex from Y marker presence
    if y_markers_found > 3:
        results["sex_determination"] = "male"
    elif y_markers_found == 0:
        results["sex_determination"] = "likely_female"
        results["haplogroup"] = "N/A (female)"
        results["confidence"] = "N/A"
        return results
    
    # Determine most likely haplogroup
    if branch_hits:
        sorted_hgs = sorted(branch_hits.items(), key=lambda x: len(x[1]), reverse=True)
        top_hg = sorted_hgs[0][0]
        
        results["haplogroup"] = top_hg
        results["branch_evidence"] = sorted_hgs[0][1]
        
        # ALWAYS low confidence for consumer arrays
        marker_count = len(sorted_hgs[0][1])
        if marker_count >= 2:
            results["confidence"] = "low - possibly correct major haplogroup"
        else:
            results["confidence"] = "very low - single marker only"
        
        if len(sorted_hgs) > 1:
            results["possible_subclades"] = [hg for hg, _ in sorted_hgs[1:4]]
        
        if top_hg in HAPLOGROUP_HISTORY:
            results["history"] = HAPLOGROUP_HISTORY[top_hg]
    
    return results


def analyze_haplogroups(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Complete haplogroup analysis including both mtDNA and Y-DNA.
    
    ⚠️ LOW CONFIDENCE - Consumer arrays cannot reliably determine haplogroups.
    These results are INDICATIVE ONLY and should not be used for:
    - Genealogical matching
    - Family tree research
    - Ancestry documentation
    
    For accurate haplogroups, use dedicated testing services.
    
    Args:
        genotypes: Dict mapping rsid -> genotype
        
    Returns:
        Dict with complete haplogroup analysis
    """
    mtdna = determine_mtdna_haplogroup(genotypes)
    ydna = determine_y_haplogroup(genotypes)
    
    return {
        "disclaimer": HAPLOGROUP_DISCLAIMER,
        "mtDNA": {
            "haplogroup": mtdna["haplogroup"],
            "confidence": mtdna["confidence"],
            "confidence_note": mtdna["confidence_note"],
            "markers_found": mtdna["markers_found"],
            "lineage": "maternal",
            "history": mtdna.get("history"),
            "evidence": mtdna.get("branch_evidence", []),
            "recommendation": mtdna["recommendation"],
            "pmid": mtdna.get("pmid", [])
        },
        "Y_DNA": {
            "haplogroup": ydna["haplogroup"],
            "confidence": ydna["confidence"],
            "confidence_note": ydna.get("confidence_note", ""),
            "markers_found": ydna["markers_found"],
            "lineage": "paternal",
            "sex": ydna["sex_determination"],
            "history": ydna.get("history"),
            "evidence": ydna.get("branch_evidence", []),
            "recommendation": ydna["recommendation"],
            "pmid": ydna.get("pmid", [])
        },
        "summary": _generate_haplogroup_summary(mtdna, ydna),
        "methodology": {
            "description": "Haplogroup inference from consumer DNA array markers",
            "limitations": [
                "Consumer arrays test only ~40 haplogroup markers total",
                "The mtDNA and Y-DNA phylogenetic trees contain thousands of defining SNPs",
                "Missing branch-defining SNPs can lead to incorrect classification",
                "Subclade resolution is not possible with array data",
                "Results may conflict with dedicated haplogroup testing"
            ],
            "recommendation": "For genealogical research or accurate haplogroup determination, "
                            "use FTDNA Big Y / mtDNA Full Sequence, YFull, or whole genome sequencing",
            "pmid": ["19165223", "18285812", "27654910"]
        }
    }


def _generate_haplogroup_summary(mtdna: Dict, ydna: Dict) -> str:
    """Generate human-readable haplogroup summary."""
    lines = []
    
    lines.append("⚠️ LOW CONFIDENCE - Consumer Array Limitation")
    lines.append("These haplogroup calls are INDICATIVE ONLY. Use dedicated testing for accuracy.")
    lines.append("")
    
    if mtdna["haplogroup"] != "Unknown":
        lines.append(f"Maternal lineage (mtDNA): Possibly Haplogroup {mtdna['haplogroup']}")
        lines.append(f"  Confidence: {mtdna['confidence']}")
        if mtdna.get("history"):
            hist = mtdna["history"]
            lines.append(f"  Origin: {hist.get('origin', 'Unknown')}")
            if hist.get("migration"):
                lines.append(f"  History: {hist['migration'][0]}")
        lines.append(f"  For accurate result: {mtdna['recommendation']}")
    else:
        lines.append("Maternal lineage: Unable to determine (insufficient markers)")
    
    lines.append("")
    
    if ydna["sex_determination"] == "likely_female":
        lines.append("Paternal lineage (Y-DNA): Not applicable (no Y-chromosome detected)")
    elif ydna["haplogroup"] != "Unknown":
        lines.append(f"Paternal lineage (Y-DNA): Possibly Haplogroup {ydna['haplogroup']}")
        lines.append(f"  Confidence: {ydna['confidence']}")
        if ydna.get("history"):
            hist = ydna["history"]
            lines.append(f"  Origin: {hist.get('origin', 'Unknown')}")
            if hist.get("migration"):
                lines.append(f"  History: {hist['migration'][0]}")
        lines.append(f"  For accurate result: {ydna['recommendation']}")
    else:
        lines.append("Paternal lineage: Unable to determine (insufficient markers)")
    
    return "\n".join(lines)
