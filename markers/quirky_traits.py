"""
Quirky/Fun Traits Genetics v5.0

The entertaining genetics:
- Earwax type / body odor (ABCC11)
- Asparagus urine smell (OR2M7)
- Photic sneeze reflex
- Cilantro soapy taste (OR6A2)
- Mosquito attraction
- Morning breath
- Motion sickness
- Fear of heights genetics
- Coffee taste
- Misophonia (sound sensitivity)
- And more!

These are for fun and self-discovery - no medical implications.
All markers with PMID references.
"""

from typing import Dict, List, Any

# =============================================================================
# EARWAX & BODY ODOR (ABCC11)
# =============================================================================

ABCC11_MARKERS = {
    "rs17822931": {
        "gene": "ABCC11",
        "variant": "Gly180Arg",
        "function": "ATP-binding cassette transporter - apocrine secretion",
        "risk_allele": "T",  # Dry earwax
        "frequency": {"EUR": 0.10, "AFR": 0.01, "EAS": 0.90, "AMR": 0.50},
        "effect": {
            "TT": "Dry earwax, REDUCED body odor - common in East Asians",
            "CC": "Wet earwax, normal body odor - common in Europeans/Africans",
            "CT": "Intermediate (usually wet)"
        },
        "category": "body",
        "evidence": "definitive",
        "pmid": ["16444273", "19710689"],
        "fun_facts": [
            "Same gene controls both earwax type AND body odor",
            "TT genotype (dry earwax) = much less body odor",
            "East Asians rarely need deodorant (genetically!)",
            "This is why deodorant isn't sold as much in East Asia",
            "Also affects colostrum (breast milk) production"
        ],
        "actionable": {
            "TT": [
                "Dry, flaky earwax (not sticky)",
                "Significantly reduced body odor",
                "May not need regular deodorant use",
                "Common if East Asian ancestry"
            ],
            "CC": [
                "Wet, sticky earwax",
                "Normal body odor (may benefit from deodorant)",
                "Common in European and African ancestry"
            ]
        }
    },
}

# =============================================================================
# ASPARAGUS URINE SMELL
# =============================================================================

ASPARAGUS_MARKERS = {
    "rs4481887": {
        "gene": "OR2M7",
        "variant": "Olfactory receptor near asparagus metabolite detection",
        "function": "Smell detection for asparagus metabolites",
        "risk_allele": "A",
        "frequency": {"EUR": 0.50},
        "effect": {
            "AA": "Can smell asparagus metabolites in urine",
            "GG": "CANNOT smell asparagus urine - anosmia to this compound",
            "AG": "Intermediate ability"
        },
        "category": "smell",
        "evidence": "strong",
        "pmid": ["27128349", "26197943"],
        "fun_facts": [
            "Almost everyone produces smelly asparagus urine",
            "But only some people can SMELL it",
            "If you don't notice it, you're likely GG genotype",
            "Not the ability to produce - the ability to detect!",
            "~40% of Europeans can't smell it"
        ],
        "actionable": {
            "GG": [
                "You probably don't notice asparagus urine smell",
                "You still produce it - just can't smell it",
                "Others might notice even if you don't!"
            ],
            "AA": [
                "You can detect asparagus urine smell",
                "Now you know why some people claim they 'don't get it'"
            ]
        }
    },
}

# =============================================================================
# PHOTIC SNEEZE REFLEX (SUN SNEEZING)
# =============================================================================

PHOTIC_SNEEZE_MARKERS = {
    "rs10427255": {
        "gene": "Near ZEB2",
        "variant": "Photic sneeze GWAS hit",
        "function": "Unknown - trigeminal nerve sensitivity?",
        "risk_allele": "C",
        "frequency": {"EUR": 0.25},
        "effect": {
            "CC": "Strong photic sneeze reflex (ACHOO syndrome)",
            "CT": "Moderate likelihood",
            "TT": "Low likelihood"
        },
        "category": "reflex",
        "evidence": "moderate",
        "pmid": ["20585627"],
        "fun_facts": [
            "Called ACHOO syndrome (Autosomal Cholinergic Helio-Ophthalmic Outburst)",
            "18-35% of population has this",
            "Triggered by looking at bright light, especially sun",
            "May involve crossed wires in trigeminal nerve",
            "Runs in families - dominant inheritance pattern"
        ],
        "actionable": {
            "CC": [
                "High probability of photic sneeze reflex",
                "You probably sneeze when going into bright sunlight",
                "Wear sunglasses to prevent",
                "Be careful when driving into sunlight!"
            ]
        }
    },
}

# =============================================================================
# CILANTRO TASTE (SOAPY)
# =============================================================================

CILANTRO_MARKERS = {
    "rs72921001": {
        "gene": "OR6A2",
        "variant": "Olfactory receptor 6A2",
        "function": "Detects aldehydes in cilantro",
        "risk_allele": "A",
        "frequency": {"EUR": 0.15, "EAS": 0.10, "SAS": 0.05, "AFR": 0.12},
        "effect": {
            "AA": "Cilantro tastes like SOAP - strong aversion",
            "GA": "May have mild soapy taste",
            "GG": "Cilantro tastes normal/pleasant"
        },
        "category": "taste",
        "evidence": "strong",
        "pmid": ["22927850"],
        "fun_facts": [
            "It's genetic - you're not just being picky!",
            "Detects aldehyde compounds that give soapy taste",
            "~14% of Europeans, ~3-7% of South/East Asians",
            "Julia Child famously hated cilantro",
            "Also called 'coriander' in many countries"
        ],
        "actionable": {
            "AA": [
                "Cilantro genuinely tastes like soap to you",
                "Not a preference - it's your olfactory genetics",
                "Use parsley, basil, or other herbs as substitute",
                "Crushing cilantro can reduce aldehydes slightly",
                "You're genetically justified in hating cilantro!"
            ],
            "GG": [
                "Cilantro tastes normal to you",
                "Now you understand why some people hate it",
                "They're not being dramatic - it's genetic"
            ]
        }
    },
}

# =============================================================================
# MOSQUITO ATTRACTION
# =============================================================================

MOSQUITO_MARKERS = {
    "rs62323883": {
        "gene": "Unknown (skin volatiles)",
        "variant": "Mosquito attraction GWAS hit",
        "function": "Affects skin compound production",
        "risk_allele": "G",
        "frequency": {"EUR": 0.30},
        "effect": {
            "GG": "More attractive to mosquitoes",
            "AA": "Less attractive"
        },
        "category": "body",
        "evidence": "moderate",
        "pmid": ["26077489", "34535578"],
        "fun_facts": [
            "Some people really are mosquito magnets",
            "CO2, body heat, lactic acid, skin microbiome all matter",
            "Blood type O may be more attractive (some evidence)",
            "Genetics explains ~67% of attractiveness variation",
            "Pregnant women more attractive (higher body temp, CO2)"
        ],
        "note": "Research is ongoing - multiple factors involved"
    },
    # Blood type (ABO) also associated
    "rs8176719": {
        "gene": "ABO",
        "variant": "Blood type O determinant",
        "function": "ABO blood group - also affects mosquito attraction",
        "risk_allele": "G",  # deletion = O
        "frequency": {"EUR": 0.40},
        "effect": {
            "note": "Type O may be more attractive to Aedes mosquitoes",
            "evidence": "Mixed - some studies support, others don't"
        },
        "category": "body",
        "evidence": "limited",
        "pmid": ["15311477"]
    },
}

# =============================================================================
# MOTION SICKNESS
# =============================================================================

MOTION_SICKNESS_MARKERS = {
    "rs8068318": {
        "gene": "GPD2",
        "variant": "Motion sickness GWAS hit",
        "function": "Glycerol-3-phosphate dehydrogenase",
        "risk_allele": "T",
        "frequency": {"EUR": 0.30},
        "effect": {
            "TT": "Higher motion sickness susceptibility"
        },
        "category": "vestibular",
        "evidence": "moderate",
        "pmid": ["25891879"],
        "fun_facts": [
            "Motion sickness is genetic - ~57% heritability",
            "More common in women",
            "Related to vestibular system sensitivity",
            "Those who get carsick often have more sensitive balance systems",
            "May be evolutionary - detecting poisoning by dizziness"
        ],
        "actionable": {
            "TT": [
                "Higher genetic susceptibility to motion sickness",
                "Sit in front seat of car",
                "Look at horizon when on boats",
                "Ginger may help",
                "Antihistamines (dramamine) effective"
            ]
        }
    },
    "rs12052276": {
        "gene": "MUTED",
        "variant": "Motion sickness variant",
        "function": "Unknown mechanism",
        "risk_allele": "C",
        "frequency": {"EUR": 0.40},
        "effect": "Associated with motion sickness",
        "category": "vestibular",
        "evidence": "moderate",
        "pmid": ["25891879"]
    },
}

# =============================================================================
# FEAR OF HEIGHTS (ACROPHOBIA)
# =============================================================================

ACROPHOBIA_MARKERS = {
    "rs3809162": {
        "gene": "TMEM132D",
        "variant": "Panic/anxiety-related variant",
        "function": "Transmembrane protein - anxiety circuits",
        "risk_allele": "T",
        "frequency": {"EUR": 0.65},
        "effect": {
            "TT": "Associated with height fear and panic disorder",
            "note": "Also associated with panic disorder in general"
        },
        "category": "psychology",
        "evidence": "moderate",
        "pmid": ["21573508", "23434960"],
        "fun_facts": [
            "Fear of heights is partially genetic (~30% heritable)",
            "Same gene variants associated with panic disorder",
            "Visual-vestibular mismatch at heights triggers fear",
            "Exposure therapy can help override genetic tendency",
            "May have been protective for ancestors"
        ],
        "actionable": {
            "TT": [
                "May be more prone to height fear/anxiety",
                "Not destiny - exposure therapy helps",
                "Cognitive behavioral therapy effective"
            ]
        }
    },
}

# =============================================================================
# OTHER FUN TRAITS
# =============================================================================

OTHER_FUN_MARKERS = {
    # Coffee consumption
    "rs4410790": {
        "gene": "AHR",
        "variant": "Coffee consumption GWAS",
        "function": "Aryl hydrocarbon receptor - caffeine signaling",
        "risk_allele": "C",
        "frequency": {"EUR": 0.45},
        "effect": {
            "CC": "Tends to drink MORE coffee",
            "TT": "Tends to drink less coffee"
        },
        "category": "consumption",
        "evidence": "strong",
        "pmid": ["21490707", "25288136"],
        "fun_facts": [
            "Coffee consumption is genetic (~45% heritability)",
            "Faster metabolizers drink more (need more for effect)",
            "Multiple genes affect coffee drinking behavior",
            "CYP1A2 affects metabolism, AHR affects consumption drive"
        ]
    },
    
    # Alcohol flush (covered in nutrition but fun)
    "rs671": {
        "gene": "ALDH2",
        "reference": "See nutrition_comprehensive.py",
        "fun_fact": "Asian flush - you either have it or you don't"
    },
    
    # Sleep duration
    "rs121912617": {
        "gene": "DEC2",
        "variant": "Short sleeper mutation (rare)",
        "function": "Circadian transcription factor",
        "risk_allele": "T",
        "frequency": {"EUR": 0.001},  # Very rare
        "effect": {
            "note": "True 'short sleepers' are very rare",
            "effect": "Can function on 4-6 hours sleep without impairment"
        },
        "category": "sleep",
        "evidence": "strong",
        "pmid": ["19679812"],
        "fun_facts": [
            "Natural short sleepers are EXTREMELY rare (<1%)",
            "Most people who claim to need little sleep are wrong",
            "They're just used to being sleep-deprived",
            "True short sleepers have this mutation",
            "Need <6h sleep with no negative effects"
        ]
    },
    
    # Unibrow
    "rs12155314": {
        "gene": "PAX3",
        "variant": "Eyebrow thickness/unibrow",
        "function": "Transcription factor - facial development",
        "risk_allele": "T",
        "frequency": {"EUR": 0.30},
        "effect": "Associated with monobrow/unibrow tendency",
        "category": "appearance",
        "evidence": "moderate",
        "pmid": ["26926045"],
        "fun_facts": [
            "Unibrow tendency is genetic",
            "Same gene affects other facial features",
            "Frida Kahlo embraced her genetic unibrow"
        ]
    },
    
    # Sneeze style
    "rs11856391": {
        "gene": "Unknown",
        "variant": "Sneeze style GWAS",
        "function": "Unknown",
        "risk_allele": "G",
        "frequency": {"EUR": 0.45},
        "effect": "Associated with sneeze volume/style",
        "category": "reflex",
        "evidence": "limited",
        "pmid": ["20585627"],
        "fun_fact": "Your sneeze style may be genetic - so is your parent's"
    },
    
    # Dimples
    "rs11850334": {
        "gene": "Near FOXL2",
        "variant": "Cheek dimple association",
        "function": "Facial muscle development",
        "risk_allele": "A",
        "frequency": {"EUR": 0.25},
        "effect": "Associated with cheek dimples",
        "category": "appearance",
        "evidence": "limited",
        "pmid": ["26926045"],
        "fun_fact": "Dimples are caused by muscle variation - genetic!"
    },
    
    # Cleft chin
    "rs1960808": {
        "gene": "Unknown",
        "variant": "Chin dimple/cleft chin",
        "function": "Facial development",
        "risk_allele": "G",
        "frequency": {"EUR": 0.30},
        "effect": "Associated with cleft/dimpled chin",
        "category": "appearance",
        "evidence": "limited",
        "pmid": ["26926045"]
    },
    
    # Toe webbing
    "rs11866098": {
        "gene": "Near ZNF804A",
        "variant": "Toe webbing tendency",
        "function": "Digit separation",
        "risk_allele": "A",
        "frequency": {"EUR": 0.15},
        "effect": "Associated with mild toe webbing (syndactyly tendency)",
        "category": "anatomy",
        "evidence": "limited",
        "pmid": ["26926045"]
    },
    
    # Perfect pitch (very polygenic)
    "rs3057": {
        "gene": "SLC6A4",
        "variant": "Associated with musical ability",
        "function": "Serotonin transporter",
        "risk_allele": "A",
        "frequency": {"EUR": 0.40},
        "effect": "Weakly associated with absolute pitch",
        "category": "ability",
        "evidence": "limited",
        "pmid": ["19664519"],
        "fun_facts": [
            "Perfect pitch is genetic but also requires early training",
            "More common in speakers of tonal languages",
            "Very polygenic - many genes contribute",
            "Only 1 in 10,000 people have true perfect pitch"
        ]
    },
}

# =============================================================================
# COMBINE ALL QUIRKY MARKERS
# =============================================================================

QUIRKY_TRAITS_MARKERS = {
    **ABCC11_MARKERS,
    **ASPARAGUS_MARKERS,
    **PHOTIC_SNEEZE_MARKERS,
    **CILANTRO_MARKERS,
    **MOSQUITO_MARKERS,
    **MOTION_SICKNESS_MARKERS,
    **ACROPHOBIA_MARKERS,
    **OTHER_FUN_MARKERS,
}

# =============================================================================
# ANALYSIS FUNCTION
# =============================================================================

def generate_quirky_traits_report(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """Generate fun quirky traits report."""
    
    results = []
    
    # Earwax/body odor
    abcc11 = genotypes.get("rs17822931", "CC")
    if abcc11 == "TT":
        results.append({
            "trait": "Earwax & Body Odor",
            "result": "Dry earwax, reduced body odor",
            "fun_fact": "You may not need deodorant - common in East Asian ancestry!"
        })
    elif abcc11 == "CC":
        results.append({
            "trait": "Earwax & Body Odor",
            "result": "Wet earwax, normal body odor",
            "fun_fact": "Common in European and African ancestry"
        })
    
    # Asparagus
    asp = genotypes.get("rs4481887", "AG")
    if asp == "GG":
        results.append({
            "trait": "Asparagus Urine Smell",
            "result": "Cannot smell it",
            "fun_fact": "You produce it, you just can't detect it!"
        })
    elif asp == "AA":
        results.append({
            "trait": "Asparagus Urine Smell",
            "result": "Can smell it",
            "fun_fact": "~60% of people can detect this distinctive odor"
        })
    
    # Photic sneeze
    photic = genotypes.get("rs10427255", "TT")
    if photic == "CC":
        results.append({
            "trait": "Photic Sneeze Reflex",
            "result": "Likely sun sneezer (ACHOO syndrome)",
            "fun_fact": "Bright light triggers your sneezes - it's genetic!"
        })
    
    # Cilantro
    cil = genotypes.get("rs72921001", "GG")
    if cil == "AA":
        results.append({
            "trait": "Cilantro Taste",
            "result": "Tastes like SOAP 🧼",
            "fun_fact": "You're genetically justified - it's not picky eating!"
        })
    elif cil == "GG":
        results.append({
            "trait": "Cilantro Taste",
            "result": "Normal/pleasant taste",
            "fun_fact": "Those who hate it aren't being dramatic - they taste it differently"
        })
    
    # Motion sickness
    motion = genotypes.get("rs8068318", "CC")
    if motion == "TT":
        results.append({
            "trait": "Motion Sickness",
            "result": "Higher susceptibility",
            "fun_fact": "Your vestibular system is genetically more sensitive"
        })
    
    return {
        "quirky_traits": results,
        "markers_analyzed": sum(1 for rs in QUIRKY_TRAITS_MARKERS if rs in genotypes),
        "note": "These are for fun and self-discovery - no medical implications!"
    }

# Export
__all__ = [
    'QUIRKY_TRAITS_MARKERS',
    'ABCC11_MARKERS',
    'ASPARAGUS_MARKERS',
    'PHOTIC_SNEEZE_MARKERS',
    'CILANTRO_MARKERS',
    'MOSQUITO_MARKERS',
    'MOTION_SICKNESS_MARKERS',
    'ACROPHOBIA_MARKERS',
    'OTHER_FUN_MARKERS',
    'generate_quirky_traits_report',
]
