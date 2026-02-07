"""
Trait Predictions with Statistical Rigor

Provides trait predictions with:
- Odds ratios with 95% confidence intervals
- Effect sizes with standard errors
- Probability estimates with uncertainty
- Confidence levels

Author: OpenClaw AI
Date: 2026-02-07
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
import sys
import math

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from personal_genomics.statistics import (
        ConfidenceLevel,
        trait_probability_ci,
        effect_size_ci,
        confidence_to_color,
    )
    STATS_AVAILABLE = True
except ImportError:
    STATS_AVAILABLE = False
    class ConfidenceLevel:
        DEFINITIVE = "DEFINITIVE"
        HIGH = "HIGH"
        MEDIUM = "MEDIUM"
        LOW = "LOW"
        UNCERTAIN = "UNCERTAIN"


# =============================================================================
# TRAIT EFFECT SIZES (from GWAS)
# =============================================================================

# Odds ratios and effect sizes from published GWAS
TRAIT_EFFECTS = {
    "rs12913832": {
        "gene": "HERC2/OCA2",
        "trait": "Blue eyes",
        "baseline_probability": 0.20,  # ~20% of world population
        "effect_allele": "G",
        "odds_ratio": 14.0,  # Very strong effect
        "odds_ratio_ci_lower": 10.5,
        "odds_ratio_ci_upper": 18.7,
        "pmid": "18172690",
        "interpretation": {
            "GG": "Blue/gray eyes very likely",
            "AG": "Green/hazel likely, blue possible",
            "AA": "Brown eyes very likely"
        }
    },
    "rs16891982": {
        "gene": "SLC45A2",
        "trait": "Light skin",
        "baseline_probability": 0.30,
        "effect_allele": "G",
        "odds_ratio": 6.0,
        "odds_ratio_ci_lower": 4.2,
        "odds_ratio_ci_upper": 8.6,
        "pmid": "16357253",
    },
    "rs1426654": {
        "gene": "SLC24A5",
        "trait": "Light skin",
        "baseline_probability": 0.30,
        "effect_allele": "A",
        "odds_ratio": 35.0,  # Very strong effect, nearly fixed in Europeans
        "odds_ratio_ci_lower": 25.0,
        "odds_ratio_ci_upper": 49.0,
        "pmid": "16357253",
    },
    "rs1805007": {
        "gene": "MC1R",
        "trait": "Red hair",
        "baseline_probability": 0.02,  # ~2% of world
        "effect_allele": "T",
        "odds_ratio": 18.0,
        "odds_ratio_ci_lower": 12.0,
        "odds_ratio_ci_upper": 27.0,
        "pmid": "11260714",
    },
    "rs4988235": {
        "gene": "LCT/MCM6",
        "trait": "Lactase persistence",
        "baseline_probability": 0.35,  # ~35% globally
        "effect_allele": "A",
        "odds_ratio": 50.0,  # Near-complete penetrance
        "odds_ratio_ci_lower": 35.0,
        "odds_ratio_ci_upper": 71.0,
        "pmid": "14681826",
    },
    "rs17822931": {
        "gene": "ABCC11",
        "trait": "Wet earwax",
        "baseline_probability": 0.70,  # Wet is ancestral/common
        "effect_allele": "C",  # C = wet
        "odds_ratio": 25.0,
        "odds_ratio_ci_lower": 18.0,
        "odds_ratio_ci_upper": 35.0,
        "pmid": "16444273",
    },
    "rs671": {
        "gene": "ALDH2",
        "trait": "Alcohol flush reaction",
        "baseline_probability": 0.08,  # ~8% globally
        "effect_allele": "A",
        "odds_ratio": 100.0,  # Near-complete penetrance
        "odds_ratio_ci_lower": 70.0,
        "odds_ratio_ci_upper": 143.0,
        "pmid": "10379521",
    },
    "rs762551": {
        "gene": "CYP1A2",
        "trait": "Slow caffeine metabolism",
        "baseline_probability": 0.30,
        "effect_allele": "C",  # C = slow
        "odds_ratio": 2.5,
        "odds_ratio_ci_lower": 1.8,
        "odds_ratio_ci_upper": 3.5,
        "pmid": "16522833",
    },
    "rs4680": {
        "gene": "COMT",
        "trait": "Worrier phenotype (Val158Met)",
        "baseline_probability": 0.25,
        "effect_allele": "A",  # Met allele
        "odds_ratio": 1.8,
        "odds_ratio_ci_lower": 1.4,
        "odds_ratio_ci_upper": 2.3,
        "pmid": "18227835",
    },
    "rs53576": {
        "gene": "OXTR",
        "trait": "Enhanced empathy/social behavior",
        "baseline_probability": 0.40,
        "effect_allele": "G",
        "odds_ratio": 1.5,
        "odds_ratio_ci_lower": 1.2,
        "odds_ratio_ci_upper": 1.9,
        "pmid": "19934046",
    },
}


@dataclass
class TraitPrediction:
    """Trait prediction with statistical confidence."""
    trait: str
    gene: str
    rsid: str
    genotype: str
    probability: float
    ci_lower: float
    ci_upper: float
    odds_ratio: float
    or_ci_lower: float
    or_ci_upper: float
    confidence: str
    effect_allele: str
    effect_count: int
    based_on: str
    interpretation: str
    pmid: str = ""
    
    def to_dict(self) -> dict:
        return {
            "trait": self.trait,
            "gene": self.gene,
            "rsid": self.rsid,
            "genotype": self.genotype,
            "probability": round(self.probability, 3),
            "ci_lower": round(self.ci_lower, 3),
            "ci_upper": round(self.ci_upper, 3),
            "odds_ratio": round(self.odds_ratio, 2),
            "or_ci_lower": round(self.or_ci_lower, 2),
            "or_ci_upper": round(self.or_ci_upper, 2),
            "confidence": self.confidence,
            "based_on": self.based_on,
            "interpretation": self.interpretation,
            "pmid": self.pmid,
        }


def predict_trait(
    rsid: str,
    genotype: str,
    effect_data: Dict = None
) -> Optional[TraitPrediction]:
    """
    Predict trait with statistical confidence interval.
    
    Args:
        rsid: SNP ID
        genotype: User's genotype (e.g., "AG")
        effect_data: Optional override for effect data
        
    Returns:
        TraitPrediction with probability and CI, or None if no data
    """
    data = effect_data or TRAIT_EFFECTS.get(rsid)
    if not data:
        return None
    
    effect_allele = data.get("effect_allele", "")
    baseline = data.get("baseline_probability", 0.5)
    odds_ratio = data.get("odds_ratio", 1.0)
    or_lower = data.get("odds_ratio_ci_lower", odds_ratio * 0.7)
    or_upper = data.get("odds_ratio_ci_upper", odds_ratio * 1.4)
    
    # Count effect alleles
    effect_count = genotype.upper().count(effect_allele.upper())
    
    # Calculate probability using odds ratio
    baseline_odds = baseline / (1 - baseline) if baseline < 1 else 100
    adjusted_odds = baseline_odds * (odds_ratio ** effect_count)
    probability = adjusted_odds / (1 + adjusted_odds)
    
    # Calculate CI bounds
    lower_odds = baseline_odds * (or_lower ** effect_count)
    upper_odds = baseline_odds * (or_upper ** effect_count)
    ci_lower = lower_odds / (1 + lower_odds)
    ci_upper = upper_odds / (1 + upper_odds)
    
    # Cap at valid range
    probability = max(0.001, min(0.999, probability))
    ci_lower = max(0.001, min(0.999, ci_lower))
    ci_upper = max(0.001, min(0.999, ci_upper))
    
    # Determine confidence level
    if or_lower > 1.0 or or_upper < 1.0:
        # CI excludes 1.0 - statistically significant
        if odds_ratio >= 5.0:
            conf_level = ConfidenceLevel.HIGH
        else:
            conf_level = ConfidenceLevel.MEDIUM
    else:
        conf_level = ConfidenceLevel.LOW
    
    # Get interpretation
    interpretations = data.get("interpretation", {})
    interpretation = interpretations.get(genotype, f"{effect_count} effect allele(s)")
    
    return TraitPrediction(
        trait=data.get("trait", "Unknown"),
        gene=data.get("gene", ""),
        rsid=rsid,
        genotype=genotype,
        probability=probability,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        odds_ratio=odds_ratio,
        or_ci_lower=or_lower,
        or_ci_upper=or_upper,
        confidence=conf_level.value if hasattr(conf_level, 'value') else str(conf_level),
        effect_allele=effect_allele,
        effect_count=effect_count,
        based_on=f"{data.get('gene', '')} {rsid} {genotype}",
        interpretation=interpretation,
        pmid=data.get("pmid", ""),
    )


def predict_all_traits(genotypes: Dict[str, str]) -> Dict[str, TraitPrediction]:
    """
    Predict all traits with statistical confidence.
    
    Args:
        genotypes: Dict of rsid -> genotype
        
    Returns:
        Dict of rsid -> TraitPrediction
    """
    results = {}
    
    for rsid in TRAIT_EFFECTS.keys():
        if rsid in genotypes:
            geno = genotypes[rsid]
            if geno and geno not in ("--", "00", "??"):
                prediction = predict_trait(rsid, geno)
                if prediction:
                    results[rsid] = prediction
    
    return results


def get_trait_summary(genotypes: Dict[str, str]) -> Dict[str, Any]:
    """
    Get summary of trait predictions for dashboard.
    
    Args:
        genotypes: Dict of rsid -> genotype
        
    Returns:
        Dict with all trait predictions formatted for display
    """
    predictions = predict_all_traits(genotypes)
    
    summary = {
        "predictions": [],
        "high_confidence": [],
        "traits_analyzed": len(predictions),
        "traits_available": len(TRAIT_EFFECTS),
    }
    
    for rsid, pred in predictions.items():
        pred_dict = pred.to_dict()
        summary["predictions"].append(pred_dict)
        
        if pred.confidence in ("HIGH", "DEFINITIVE"):
            summary["high_confidence"].append(pred_dict)
    
    # Sort by confidence and probability
    summary["predictions"].sort(
        key=lambda x: (x["confidence"] != "HIGH", -x["probability"]),
    )
    
    return summary
