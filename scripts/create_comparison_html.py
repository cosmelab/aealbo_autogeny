#!/usr/bin/env python3
"""
Create embedded HTML comparison report with actual figures
"""
import base64
import os
from pathlib import Path

PROJECT_DIR = Path("$HOME/projects/albopictus-autogeny")

def encode_image(path):
    """Encode image to base64"""
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode()
    ext = path.split(".")[-1].lower()
    if ext == "jpg":
        ext = "jpeg"
    return f"data:image/{ext};base64,{data}"

# Define figure pairs (original, replicated)
figures = {
    "Figure 1a - 4-way Venn (all comparisons)": {
        "manuscript": str(PROJECT_DIR / "manuscript/marker_output/FIgures/_page_0_Figure_0.jpeg"),
        "replicated": str(PROJECT_DIR / "output/figures_comparison/replicated/4_way_venn_significant_snps.png"),
        "description": "4-way Venn: NEW_MAN_AUT, MAN_AUT, NEW_AUT, NEW_MAN",
        "status": "ORIGINAL - NEEDS CONTAINER",
        "notes": "17 SNPs common across all comparisons (ggvenn package missing)"
    },
    "Figure 1b - 2-way Venn (3-pop PCAdapt vs Outflank)": {
        "replicated": str(PROJECT_DIR / "output/figures_comparison/replicated/venn_3_pop_pcadapt_outflank.png"),
        "description": "PCAdapt vs Outflank for 3-population analysis",
        "status": "ORIGINAL - NEEDS CONTAINER",
        "notes": "SNPs identified by both methods in 3-pop comparison"
    },
    "Figure 1c - 2-way Venn (NEW vs AUT)": {
        "replicated": str(PROJECT_DIR / "output/figures_comparison/replicated/venn_NEW_AUT_pcadapt_outflank.png"),
        "description": "PCAdapt vs Outflank for NEW vs AUT comparison",
        "status": "ORIGINAL - NEEDS CONTAINER",
        "notes": "Non-autogenous field vs Autogenous"
    },
    "Figure 1d - 2-way Venn (MAN vs AUT)": {
        "replicated": str(PROJECT_DIR / "output/figures_comparison/replicated/venn_MAN_AUT_pcadapt_outflank.png"),
        "description": "PCAdapt vs Outflank for MAN vs AUT comparison",
        "status": "ORIGINAL - NEEDS CONTAINER",
        "notes": "Non-autogenous lab vs Autogenous"
    },
    "Figure 2 - Genotype Frequencies": {
        "manuscript": str(PROJECT_DIR / "manuscript/marker_output/FIgures/_page_2_Figure_0.jpeg"),
        "replicated": str(PROJECT_DIR / "output/figures_comparison/replicated/significant_17_snps_genotypes.png"),
        "description": "Genotype frequencies for 17 outlier SNPs",
        "status": "VALIDATED - EXACT MATCH",
        "notes": "File sizes identical: 10,774 bytes"
    },
    "Figure 3 - Cluster Bar Plots": {
        "manuscript": str(PROJECT_DIR / "manuscript/marker_output/FIgures/_page_4_Figure_0.jpeg"),
        "replicated": str(PROJECT_DIR / "output/figures_comparison/replicated/figure_3_cluster_bars.png"),
        "description": "Linkage cluster sizes for AUTO and NON-AUTO-FIELD populations",
        "status": "GENERATED",
        "notes": "14 clusters across 3 chromosomes (from Rmd lines 4205-4245)"
    },
    "Figure 4 - Cluster 14 Scaffold": {
        "manuscript": str(PROJECT_DIR / "manuscript/marker_output/FIgures/_page_6_Figure_0.jpeg"),
        "original": str(PROJECT_DIR / "dropbox_original/autogenous/output/ldna/figures/cluster_14_aut_scaffolds.pdf"),
        "replicated": str(PROJECT_DIR / "output/ldna/figures_matched/cluster_14_aut_scaffolds.png"),
        "description": "Scaffolds containing cluster 14 SNPs",
        "status": "GENERATED",
        "notes": "Correct scaffolds: 2.22, 2.27, 2.32, 2.36, 2.38, 2.40"
    },
    "Figure 5 - Fst Sliding Windows": {
        "manuscript": str(PROJECT_DIR / "manuscript/marker_output/FIgures/_page_8_Figure_0.jpeg"),
        "replicated": str(PROJECT_DIR / "output/fst/figures/pairwise_estimates_windows.png"),
        "description": "Pairwise Fst between AUTO and wild populations",
        "status": "VALIDATED",
        "notes": "Data identical between OLD and NEW"
    },
    "Figure 6 - Sex Fst": {
        "manuscript": str(PROJECT_DIR / "manuscript/marker_output/FIgures/_page_10_Figure_0.jpeg"),
        "replicated": str(PROJECT_DIR / "output/figures_comparison/replicated/sex_pairwise_estimates_windows.png"),
        "description": "Sex-specific Fst between males and females",
        "status": "GENERATED",
        "notes": "Sliding window Fst for sex-linked regions"
    }
}

# Build HTML
html = """<!DOCTYPE html>
<html>
<head>
    <title>Figure Comparison Report - Albopictus Autogeny GWAS</title>
    <style>
        body { font-family: 'Segoe UI', Arial, sans-serif; margin: 0; padding: 20px; background: #f0f2f5; }
        .container { max-width: 1600px; margin: 0 auto; }
        h1 { color: #1a1a2e; border-bottom: 3px solid #4CAF50; padding-bottom: 15px; }
        h2 { color: #16213e; margin-top: 40px; background: #e8f5e9; padding: 15px; border-radius: 8px; }
        .comparison { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin: 20px 0; }
        .figure-box { background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
        .figure-box h3 { margin-top: 0; color: #333; font-size: 16px; border-bottom: 2px solid #eee; padding-bottom: 10px; }
        .figure-box img { max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }
        .status { padding: 8px 16px; border-radius: 20px; display: inline-block; margin: 10px 0; font-weight: bold; }
        .validated { background: #4CAF50; color: white; }
        .generated { background: #2196F3; color: white; }
        .pending { background: #FF9800; color: white; }
        .original { background: #9C27B0; color: white; }
        .exact-match { background: #8BC34A; color: white; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; background: white; border-radius: 8px; overflow: hidden; }
        th, td { padding: 15px; text-align: left; border-bottom: 1px solid #eee; }
        th { background: #4CAF50; color: white; }
        tr:hover { background: #f5f5f5; }
        .notes { font-size: 14px; color: #666; background: #fff3e0; padding: 10px; border-radius: 4px; margin-top: 10px; }
        .description { font-size: 14px; color: #555; margin-bottom: 15px; }
        .full-width { grid-column: 1 / -1; }
    </style>
</head>
<body>
<div class="container">
    <h1>Manuscript Figure Comparison Report</h1>
    <p><strong>Project:</strong> Albopictus Autogeny GWAS - Rapid Evolution of Autogeny</p>
    <p><strong>Generated:</strong> 2026-01-06 (Labels Updated)</p>
    <p><strong>Purpose:</strong> Compare manuscript figures with replicated outputs</p>

    <h2>Summary Table</h2>
    <table>
        <tr>
            <th>Figure</th>
            <th>Description</th>
            <th>Status</th>
            <th>Notes</th>
        </tr>
"""

# Add summary table rows
for name, info in figures.items():
    status_class = "validated" if "VALIDATED" in info["status"] else ("pending" if "PENDING" in info["status"] else ("original" if "ORIGINAL" in info["status"] else "generated"))
    html += f"""        <tr>
            <td><strong>{name}</strong></td>
            <td>{info['description']}</td>
            <td><span class="status {status_class}">{info['status']}</span></td>
            <td>{info['notes']}</td>
        </tr>
"""

html += """    </table>
"""

# Add individual figure comparisons
for name, info in figures.items():
    status_class = "validated" if "VALIDATED" in info["status"] else ("pending" if "PENDING" in info["status"] else ("original" if "ORIGINAL" in info["status"] else "generated"))

    html += f"""
    <h2>{name}</h2>
    <p class="description">{info['description']}</p>
    <p><span class="status {status_class}">{info['status']}</span></p>
    <div class="notes">{info['notes']}</div>
    <div class="comparison">
"""

    # Manuscript version
    if "manuscript" in info and os.path.exists(info["manuscript"]):
        b64 = encode_image(info["manuscript"])
        if b64:
            html += f"""        <div class="figure-box">
            <h3>MANUSCRIPT VERSION</h3>
            <img src="{b64}" alt="{name} - Manuscript">
        </div>
"""

    # Replicated version
    if "replicated" in info and os.path.exists(info["replicated"]):
        b64 = encode_image(info["replicated"])
        if b64:
            html += f"""        <div class="figure-box">
            <h3>REPLICATED VERSION</h3>
            <img src="{b64}" alt="{name} - Replicated">
        </div>
"""
    elif "replicated_pdf" in info:
        html += f"""        <div class="figure-box">
            <h3>REPLICATED VERSION</h3>
            <p>PDF file: {os.path.basename(info['replicated_pdf'])}</p>
            <p><em>PDF cannot be embedded - view file directly</em></p>
        </div>
"""
    else:
        html += """        <div class="figure-box">
            <h3>REPLICATED VERSION</h3>
            <p><em>Not yet generated</em></p>
        </div>
"""

    html += """    </div>
"""

# Add key findings section
html += """
    <h2>Key Findings</h2>
    <table>
        <tr>
            <th>Finding</th>
            <th>Value</th>
            <th>Significance</th>
        </tr>
        <tr>
            <td>Cluster 14</td>
            <td>~7Mb on Chromosome 2</td>
            <td>Major autogeny-associated region</td>
        </tr>
        <tr>
            <td>Cluster 14 SNPs</td>
            <td>316 SNPs</td>
            <td>Distributed across 6 scaffolds</td>
        </tr>
        <tr>
            <td>Cluster 14 Scaffolds</td>
            <td>2.22, 2.27, 2.32, 2.36, 2.38, 2.40</td>
            <td>All correctly identified in replication</td>
        </tr>
        <tr>
            <td>Outlier SNPs</td>
            <td>157 total</td>
            <td>17 common across all comparisons</td>
        </tr>
        <tr>
            <td>Fst Peak</td>
            <td>Chromosome 2</td>
            <td>Co-localizes with cluster 14</td>
        </tr>
    </table>

    <h2>LDna Concepts</h2>
    <table>
        <tr>
            <th>Term</th>
            <th>Definition</th>
        </tr>
        <tr>
            <td><strong>SOC</strong></td>
            <td>Single Outlier Cluster - isolated high-LD region at branch tip</td>
        </tr>
        <tr>
            <td><strong>COC</strong></td>
            <td>Compound Outlier Cluster - complex structure with nested clusters</td>
        </tr>
        <tr>
            <td><strong>extractBranches()</strong></td>
            <td>Current LDna v2.0+ function (replaces deprecated extractClusters)</td>
        </tr>
    </table>

    <h2>Working Scripts</h2>
    <table>
        <tr>
            <th>Figure</th>
            <th>Script</th>
        </tr>
        <tr>
            <td>Figure 1 (Venn)</td>
            <td>scripts/12_ldna_figures/02_venn_diagrams.R</td>
        </tr>
        <tr>
            <td>Figure 2 (Genotypes)</td>
            <td>scripts/12_ldna_figures/03_genotype_frequencies.R</td>
        </tr>
        <tr>
            <td>Figure 4 (Scaffolds)</td>
            <td>scripts/08_figures/ldna_scaffold_plots_v2.R</td>
        </tr>
        <tr>
            <td>Figure 5/6 (Fst)</td>
            <td>scripts/04_diversity/sliding_window_fst.sh</td>
        </tr>
    </table>
</div>
</body>
</html>
"""

# Write HTML file
output_path = PROJECT_DIR / "output/manuscript_figure_comparison.html"
with open(output_path, "w") as f:
    f.write(html)

print(f"Created: {output_path}")
print(f"File size: {os.path.getsize(output_path):,} bytes")
