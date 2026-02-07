#!/usr/bin/env python3
"""
Generate comprehensive HTML dashboard for DNA analysis results.

Author: OpenClaw AI
Date: 2026-02-07
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict


def generate_dashboard(results_file: str, output_file: str) -> None:
    """Generate HTML dashboard from analysis results."""
    
    with open(results_file) as f:
        results = json.load(f)
    
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Kyle's Comprehensive DNA Analysis - {datetime.now().strftime("%B %d, %Y")}</title>
    <style>
        :root {{
            --primary: #2563eb;
            --success: #16a34a;
            --warning: #d97706;
            --danger: #dc2626;
            --bg: #f8fafc;
            --card-bg: #ffffff;
            --text: #1e293b;
            --text-muted: #64748b;
            --border: #e2e8f0;
        }}
        
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
            padding: 2rem;
        }}
        
        .container {{ max-width: 1400px; margin: 0 auto; }}
        
        header {{
            text-align: center;
            padding: 2rem 0;
            margin-bottom: 2rem;
            border-bottom: 2px solid var(--border);
        }}
        
        header h1 {{
            font-size: 2.5rem;
            background: linear-gradient(135deg, var(--primary), #7c3aed);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }}
        
        .meta-info {{
            display: flex;
            justify-content: center;
            gap: 2rem;
            margin-top: 1rem;
            color: var(--text-muted);
            font-size: 0.9rem;
        }}
        
        .section {{
            background: var(--card-bg);
            border-radius: 12px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        
        .section h2 {{
            color: var(--primary);
            margin-bottom: 1rem;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid var(--border);
        }}
        
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 1.5rem;
        }}
        
        .card {{
            background: var(--bg);
            border-radius: 8px;
            padding: 1rem;
        }}
        
        .card h3 {{
            font-size: 1rem;
            color: var(--text-muted);
            margin-bottom: 0.5rem;
        }}
        
        .bar-chart {{
            margin-top: 1rem;
        }}
        
        .bar-item {{
            display: flex;
            align-items: center;
            margin-bottom: 0.5rem;
        }}
        
        .bar-label {{
            width: 150px;
            font-size: 0.85rem;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
        
        .bar-container {{
            flex: 1;
            height: 20px;
            background: var(--border);
            border-radius: 4px;
            overflow: hidden;
            margin: 0 0.5rem;
        }}
        
        .bar {{
            height: 100%;
            background: linear-gradient(90deg, var(--primary), #7c3aed);
            border-radius: 4px;
            transition: width 0.5s ease;
        }}
        
        .bar-value {{
            width: 50px;
            text-align: right;
            font-size: 0.85rem;
            font-weight: 600;
        }}
        
        .signal {{
            display: inline-flex;
            align-items: center;
            gap: 0.5rem;
            padding: 0.5rem 1rem;
            border-radius: 20px;
            font-size: 0.85rem;
            margin: 0.25rem;
        }}
        
        .signal-strong {{ background: #dcfce7; color: #166534; }}
        .signal-moderate {{ background: #fef3c7; color: #92400e; }}
        .signal-weak {{ background: #fee2e2; color: #991b1b; }}
        .signal-none {{ background: #f1f5f9; color: #64748b; }}
        
        .finding {{
            padding: 1rem;
            border-left: 4px solid var(--primary);
            background: var(--bg);
            margin-bottom: 0.5rem;
            border-radius: 0 8px 8px 0;
        }}
        
        .finding.pathogenic {{ border-color: var(--danger); }}
        .finding.actionable {{ border-color: var(--warning); }}
        
        .finding-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        
        .badge {{
            padding: 0.25rem 0.75rem;
            border-radius: 20px;
            font-size: 0.75rem;
            font-weight: 600;
            text-transform: uppercase;
        }}
        
        .badge-danger {{ background: #fee2e2; color: var(--danger); }}
        .badge-warning {{ background: #fef3c7; color: var(--warning); }}
        .badge-success {{ background: #dcfce7; color: var(--success); }}
        .badge-neutral {{ background: #f1f5f9; color: var(--text-muted); }}
        
        .prs-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
        }}
        
        .prs-card {{
            background: var(--bg);
            border-radius: 8px;
            padding: 1rem;
            text-align: center;
        }}
        
        .prs-card h4 {{ font-size: 0.85rem; color: var(--text-muted); }}
        .prs-card .percentile {{ font-size: 2rem; font-weight: 700; }}
        .prs-card .category {{ font-size: 0.85rem; }}
        
        /* Statistical confidence styling */
        .confidence-bar {{
            height: 6px;
            background: var(--border);
            border-radius: 3px;
            margin-top: 0.5rem;
            position: relative;
            overflow: visible;
        }}
        
        .confidence-range {{
            position: absolute;
            height: 100%;
            border-radius: 3px;
            transition: all 0.3s ease;
        }}
        
        .confidence-point {{
            position: absolute;
            width: 8px;
            height: 8px;
            border-radius: 50%;
            background: #1e293b;
            top: -1px;
            transform: translateX(-50%);
        }}
        
        .confidence-badge {{
            display: inline-flex;
            align-items: center;
            gap: 0.25rem;
            padding: 0.2rem 0.5rem;
            border-radius: 4px;
            font-size: 0.7rem;
            font-weight: 600;
            text-transform: uppercase;
        }}
        
        .confidence-DEFINITIVE {{ background: #dcfce7; color: #166534; }}
        .confidence-HIGH {{ background: #d1fae5; color: #065f46; }}
        .confidence-MEDIUM {{ background: #fef3c7; color: #92400e; }}
        .confidence-LOW {{ background: #fed7aa; color: #9a3412; }}
        .confidence-UNCERTAIN {{ background: #fee2e2; color: #991b1b; }}
        
        .ci-text {{
            font-size: 0.75rem;
            color: var(--text-muted);
            margin-top: 0.25rem;
        }}
        
        .stat-tooltip {{
            position: relative;
            cursor: help;
        }}
        
        .stat-tooltip:hover::after {{
            content: attr(data-tooltip);
            position: absolute;
            bottom: 100%;
            left: 50%;
            transform: translateX(-50%);
            background: #1e293b;
            color: white;
            padding: 0.5rem;
            border-radius: 4px;
            font-size: 0.75rem;
            white-space: nowrap;
            z-index: 100;
        }}
        
        .prs-ci {{
            font-size: 0.75rem;
            color: var(--text-muted);
            margin-top: 0.25rem;
        }}
        
        .trait-tags {{
            display: flex;
            flex-wrap: wrap;
            gap: 0.5rem;
        }}
        
        .trait-tag {{
            padding: 0.5rem 1rem;
            background: var(--bg);
            border-radius: 20px;
            font-size: 0.85rem;
        }}
        
        footer {{
            text-align: center;
            margin-top: 3rem;
            padding-top: 2rem;
            border-top: 2px solid var(--border);
            color: var(--text-muted);
            font-size: 0.85rem;
        }}
        
        .disclaimer {{
            background: #fef3c7;
            border: 1px solid #fcd34d;
            border-radius: 8px;
            padding: 1rem;
            margin-bottom: 1.5rem;
            font-size: 0.9rem;
        }}
        
        @media (max-width: 768px) {{
            body {{ padding: 1rem; }}
            header h1 {{ font-size: 1.75rem; }}
            .meta-info {{ flex-direction: column; gap: 0.5rem; }}
            .bar-label {{ width: 100px; font-size: 0.75rem; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧬 Comprehensive DNA Analysis</h1>
            <div class="meta-info">
                <span>📅 {datetime.fromisoformat(results["analysis_date"]).strftime("%B %d, %Y")}</span>
                <span>📊 {results["total_snps"]:,} SNPs Analyzed</span>
                <span>🔬 9 Reference Datasets</span>
            </div>
        </header>
        
        <div class="disclaimer">
            <strong>⚠️ Important:</strong> This analysis is for educational and genealogical purposes only. 
            It is NOT intended for medical diagnosis. Consult a healthcare provider or genetic counselor 
            for medical interpretation of genetic data.
        </div>
        
        <!-- Population Ancestry Section -->
        <section class="section">
            <h2>🌍 Population Ancestry</h2>
            <div class="grid">
                <div class="card">
                    <h3>1000 Genomes - Continental Ancestry</h3>
                    <div class="bar-chart">
'''
    
    # Add superpopulation breakdown
    for superpop, pct in sorted(results["superpopulation_breakdown"].items(), 
                                 key=lambda x: x[1], reverse=True):
        html += f'''
                        <div class="bar-item">
                            <span class="bar-label">{superpop}</span>
                            <div class="bar-container">
                                <div class="bar" style="width: {min(pct, 100)}%"></div>
                            </div>
                            <span class="bar-value">{pct:.1f}%</span>
                        </div>'''
    
    html += '''
                    </div>
                </div>
                
                <div class="card">
                    <h3>Top Similar Populations (1000 Genomes)</h3>
                    <div class="bar-chart">
'''
    
    # Add top populations
    for pop in results["top_1kg_populations"][:5]:
        html += f'''
                        <div class="bar-item">
                            <span class="bar-label" title="{pop["name"]}">{pop["population"]}</span>
                            <div class="bar-container">
                                <div class="bar" style="width: {min(pop["similarity"], 100)}%"></div>
                            </div>
                            <span class="bar-value">{pop["similarity"]:.1f}%</span>
                        </div>'''
    
    html += '''
                    </div>
                </div>
                
                <div class="card">
                    <h3>Top Similar Populations (HGDP)</h3>
                    <div class="bar-chart">
'''
    
    # Add HGDP populations
    for pop in results["top_hgdp_populations"][:5]:
        html += f'''
                        <div class="bar-item">
                            <span class="bar-label">{pop["name"]}</span>
                            <div class="bar-container">
                                <div class="bar" style="width: {min(pop["similarity"], 100)}%"></div>
                            </div>
                            <span class="bar-value">{pop["similarity"]:.1f}%</span>
                        </div>'''
    
    html += '''
                    </div>
                </div>
            </div>
        </section>
        
        <!-- Ancient Ancestry Section -->
        <section class="section">
            <h2>🏛️ Ancient Ancestry Signals</h2>
            <p style="color: var(--text-muted); margin-bottom: 1rem;">
                Detection of ancestral population signals based on curated markers from ancient DNA studies.
            </p>
            <div style="display: flex; flex-wrap: wrap; gap: 0.5rem;">
'''
    
    # Add ancient signals
    signal_class_map = {
        "STRONG": "signal-strong",
        "MODERATE": "signal-moderate", 
        "WEAK": "signal-weak",
        "NONE": "signal-none"
    }
    
    for signal in results["ancient_signals"]:
        strength = signal["signal_strength"].upper()
        css_class = signal_class_map.get(strength, "signal-none")
        emoji = "🟢" if strength == "STRONG" else "🟡" if strength == "MODERATE" else "🟠" if strength == "WEAK" else "⚪"
        html += f'''
                <div class="signal {css_class}">
                    {emoji} <strong>{signal["population"]}</strong>: {strength}
                    ({signal["markers_matched"]}/{signal["total_markers"]} markers)
                </div>'''
    
    html += f'''
            </div>
            <div class="card" style="margin-top: 1rem;">
                <h3>Neanderthal Introgression</h3>
                <p>Markers detected: <strong>{results["neanderthal_markers"]}/{results["neanderthal_total"]}</strong></p>
                <p style="font-size: 0.85rem; color: var(--text-muted);">
                    Most non-African humans carry 1-4% Neanderthal DNA from admixture ~50,000 years ago.
                </p>
            </div>
        </section>
        
        <!-- Clinical Findings Section -->
        <section class="section">
            <h2>🏥 Clinical Findings</h2>
            <div class="grid">
                <div>
                    <h3 style="color: var(--danger);">Pathogenic Variants: {results["pathogenic_count"]}</h3>
'''
    
    # Add clinical findings
    for finding in results["clinvar_findings"][:5]:
        significance = finding["significance"].lower()
        badge_class = "badge-danger" if "pathogenic" in significance else "badge-warning" if "uncertain" in significance else "badge-neutral"
        html += f'''
                    <div class="finding pathogenic">
                        <div class="finding-header">
                            <strong>{finding["rsid"]} ({finding["gene"]})</strong>
                            <span class="badge {badge_class}">{finding["significance"]}</span>
                        </div>
                        <p style="font-size: 0.85rem; color: var(--text-muted);">
                            {', '.join(finding["conditions"][:2]) if finding["conditions"] else 'Condition not specified'}
                        </p>
                    </div>'''
    
    html += '''
                </div>
            </div>
        </section>
        
        <!-- Pharmacogenomics Section -->
        <section class="section">
            <h2>💊 Pharmacogenomics</h2>
            <p style="color: var(--text-muted); margin-bottom: 1rem;">
                Genetic variants affecting drug metabolism and response. Actionable findings: <strong>{}</strong>
            </p>
            <div class="grid">
'''.format(results["actionable_pgx"])
    
    # Add pharmacogenomics
    for pgx in results["pharmacogenomics"]:
        phenotype = pgx["phenotype"]
        if "poor" in phenotype.lower() or "ultra" in phenotype.lower():
            badge_class = "badge-warning"
        elif "intermediate" in phenotype.lower():
            badge_class = "badge-neutral"
        else:
            badge_class = "badge-success"
        
        html += f'''
                <div class="card">
                    <div class="finding-header">
                        <h3>{pgx["gene"]}</h3>
                        <span class="badge {badge_class}">{phenotype}</span>
                    </div>
                    <p style="font-size: 0.85rem; margin-top: 0.5rem;">
                        Activity Score: {pgx["activity_score"]}
                    </p>
                    <p style="font-size: 0.85rem; color: var(--text-muted);">
                        Affects: {', '.join(pgx["drugs_affected"][:3])}
                    </p>
                </div>'''
    
    html += '''
            </div>
        </section>
        
        <!-- Polygenic Risk Scores Section -->
        <section class="section">
            <h2>📊 Polygenic Risk Scores</h2>
            <p style="color: var(--text-muted); margin-bottom: 1rem;">
                Combined genetic risk from multiple variants. Percentile indicates position relative to population.
            </p>
            <div class="prs-grid">
'''
    
    # Add PRS results
    for prs in results["prs_results"]:
        category = prs["risk_category"]
        if category == "High":
            color = "var(--danger)"
        elif category == "Elevated":
            color = "var(--warning)"
        else:
            color = "var(--success)"
        
        html += f'''
                <div class="prs-card">
                    <h4>{prs["trait"]}</h4>
                    <div class="percentile" style="color: {color};">{prs["percentile"]:.0f}</div>
                    <div class="category" style="color: {color};">percentile • {category}</div>
                </div>'''
    
    html += '''
            </div>
        </section>
        
        <!-- GWAS Trait Associations Section -->
        <section class="section">
            <h2>🔬 Trait Associations (GWAS)</h2>
            <p style="color: var(--text-muted); margin-bottom: 1rem;">
                Genetic variants associated with various traits from genome-wide association studies.
            </p>
            <div class="trait-tags">
'''
    
    # Add notable traits
    for trait in results["notable_traits"][:10]:
        html += f'''
                <span class="trait-tag">{trait}</span>'''
    
    html += f'''
            </div>
        </section>
        
        <footer>
            <p>Generated by Personal Genomics Skill v4.4.0</p>
            <p>Analysis performed locally - no data sent to external servers</p>
            <p>Report generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        </footer>
    </div>
</body>
</html>
'''
    
    with open(output_file, 'w') as f:
        f.write(html)
    
    print(f"Dashboard generated: {output_file}")


if __name__ == "__main__":
    results_file = Path.home() / "dna-analysis" / "reports" / "overnight_2026-02-07" / "comprehensive_analysis.json"
    output_file = Path.home() / "dna-analysis" / "reports" / "overnight_2026-02-07" / "dashboard.html"
    
    generate_dashboard(str(results_file), str(output_file))
