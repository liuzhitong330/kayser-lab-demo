"""
make_figures.py — Kayser Lab Demo
Figure 1: Sleep GWAS landscape — loci per trait (GWAS Catalog data)
Figure 2: Kayser lab cross-species validation pipeline — human gene → Drosophila ortholog → phenotype
"""

import csv, os, math

OUT = os.path.dirname(__file__)

traits = []
with open(os.path.join(OUT, "sleep_gwas_traits.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        traits.append((row["trait"], int(row["gwas_loci"]), int(row["n_studies"])))

genes = []
with open(os.path.join(OUT, "validated_genes.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        genes.append(row)

# ── Figure 1: Sleep GWAS landscape ───────────────────────────────────────────
FW, FH = 680, 400
PAD_L = 220
PAD_R = 100
PAD_T = 75
PAD_B = 60
AW = FW - PAD_L - PAD_R
AH = FH - PAD_T - PAD_B

TRAIT_COLOR = {
    "insomnia":                   "#c0392b",
    "chronotype":                 "#1a5c8a",
    "sleep duration":             "#27ae60",
    "narcolepsy":                 "#8e44ad",
    "excessive daytime sleepiness": "#e67e22",
    "sleep quality":              "#16a085",
    "hypersomnia":                "#e67e22",
    "sleep latency":              "#7f8c8d",
}

max_loci = max(t[1] for t in traits)
n = len(traits)
bar_h = min(32, (AH - (n-1)*8) // n)
gap = 8

bars1 = ""
for i, (trait, loci, studies) in enumerate(traits):
    y = PAD_T + i * (bar_h + gap)
    bw = loci / max_loci * AW
    col = TRAIT_COLOR.get(trait, "#95a5a6")

    bars1 += (f'<rect x="{PAD_L}" y="{y}" width="{bw:.1f}" height="{bar_h}" '
              f'fill="{col}" opacity="0.85" rx="2"/>')
    # Trait label
    bars1 += (f'<text x="{PAD_L - 10}" y="{y + bar_h/2 + 4:.1f}" '
              f'text-anchor="end" font-size="11" fill="#333" font-weight="600">'
              f'{trait.title()}</text>')
    # Loci count right of bar
    bars1 += (f'<text x="{PAD_L + bw + 7:.1f}" y="{y + bar_h/2 + 4:.1f}" '
              f'font-size="10" fill="{col}" font-weight="700">{loci:,}</text>')
    # Study count far right
    bars1 += (f'<text x="{FW - PAD_R + 10}" y="{y + bar_h/2 + 4:.1f}" '
              f'font-size="9" fill="#999">{studies} studies</text>')

# Axis
xt1 = ""
for v in [0, 500, 1000, 2000, 3000, 4000]:
    if v > max_loci * 1.05: break
    tx = PAD_L + v / max_loci * AW
    xt1 += (f'<line x1="{tx:.1f}" y1="{PAD_T}" x2="{tx:.1f}" y2="{PAD_T+AH}" '
            f'stroke="#eee" stroke-width="1"/>'
            f'<text x="{tx:.1f}" y="{PAD_T+AH+16}" text-anchor="middle" '
            f'font-size="9" fill="#888">{v:,}</text>')
ax1 = (f'<line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T+AH}" '
       f'stroke="#ccc" stroke-width="1"/>')

svg1 = f"""<svg viewBox="0 0 {FW} {FH}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Sleep GWAS Landscape: Known Genetic Loci by Trait
  </text>
  <text x="{FW//2}" y="40" text-anchor="middle" font-size="10" fill="#666">
    GWAS Catalog (EBI) · all genome-wide significant associations per sleep trait · July 2026
  </text>
  <text x="{FW//2}" y="56" text-anchor="middle" font-size="10" fill="#444">
    Insomnia dominates with 4,422 loci — yet most remain functionally uncharacterized
  </text>
  <text x="{PAD_L + AW/2:.0f}" y="{PAD_T+AH+34}" text-anchor="middle"
        font-size="10" fill="#555">Number of GWAS loci (genome-wide significant)</text>
  {ax1}{xt1}{bars1}
</svg>"""

with open(os.path.join(OUT, "sleep_gwas.svg"), "w") as f:
    f.write(svg1)
print("Wrote sleep_gwas.svg")


# ── Figure 2: Cross-species validation pipeline ───────────────────────────────
FW2, FH2 = 720, 400
PAD_T2, PAD_B2 = 75, 30

GENE_COLORS = ["#c0392b", "#1a5c8a", "#27ae60", "#8e44ad"]
TRAIT_COLORS = {
    "excessive daytime sleepiness": "#e67e22",
    "hypersomnia (sub-threshold)":  "#e67e22",
    "sleep disruption (TDP-43)":    "#7f8c8d",
    "insomnia":                     "#c0392b",
}

n_genes = len(genes)
row_h = (FH2 - PAD_T2 - PAD_B2) / n_genes
COL_HUMAN   = 80
COL_DROSO   = 280
COL_PHENO   = 440
COL_MECH    = 600
COL_W       = 130

rows2 = ""

# Column headers
for cx, label in [(COL_HUMAN, "Human Gene"), (COL_DROSO, "Drosophila\nOrtholog"),
                  (COL_PHENO, "Fly Phenotype"), (COL_MECH, "Mechanism")]:
    rows2 += (f'<text x="{cx + COL_W/2:.0f}" y="30" text-anchor="middle" '
              f'font-size="10" fill="#555" font-weight="700">{label.replace(chr(10)," ")}</text>')
    rows2 += (f'<line x1="{cx}" y1="36" x2="{cx+COL_W}" y2="36" '
              f'stroke="#ddd" stroke-width="1"/>')
# GWAS trait column
rows2 += (f'<text x="{COL_HUMAN + COL_W/2:.0f}" y="44" text-anchor="middle" '
          f'font-size="9" fill="#888">GWAS trait</text>')

for i, g in enumerate(genes):
    col = GENE_COLORS[i]
    y_center = PAD_T2 + i * row_h + row_h / 2
    y_top    = PAD_T2 + i * row_h + 6
    y_bot    = PAD_T2 + (i+1) * row_h - 6

    # Alternating row background
    if i % 2 == 0:
        rows2 += (f'<rect x="0" y="{PAD_T2 + i*row_h:.1f}" width="{FW2}" '
                  f'height="{row_h:.1f}" fill="#f8f8f8" rx="0"/>')

    # Human gene box
    rows2 += (f'<rect x="{COL_HUMAN}" y="{y_top:.1f}" width="{COL_W}" '
              f'height="{row_h-12:.1f}" rx="4" fill="{col}" opacity="0.12"/>')
    rows2 += (f'<text x="{COL_HUMAN + COL_W/2:.0f}" y="{y_center - 5:.1f}" '
              f'text-anchor="middle" font-size="13" font-weight="700" fill="{col}">'
              f'{g["human_gene"]}</text>')
    # GWAS trait under gene name
    trait_col = TRAIT_COLORS.get(g["sleep_trait"], "#555")
    rows2 += (f'<text x="{COL_HUMAN + COL_W/2:.0f}" y="{y_center + 11:.1f}" '
              f'text-anchor="middle" font-size="8.5" fill="{trait_col}">'
              f'{g["sleep_trait"][:22]}</text>')

    # Arrow
    arrow_x1 = COL_HUMAN + COL_W + 4
    arrow_x2 = COL_DROSO - 4
    arrow_mx  = (arrow_x1 + arrow_x2) / 2
    rows2 += (f'<line x1="{arrow_x1}" y1="{y_center:.1f}" x2="{arrow_x2}" y2="{y_center:.1f}" '
              f'stroke="{col}" stroke-width="1.5" marker-end="url(#arr{i})"/>')
    rows2 += (f'<text x="{arrow_mx:.0f}" y="{y_center - 5:.1f}" text-anchor="middle" '
              f'font-size="7.5" fill="{col}">ortholog</text>')

    # Drosophila gene box
    rows2 += (f'<rect x="{COL_DROSO}" y="{y_top:.1f}" width="{COL_W}" '
              f'height="{row_h-12:.1f}" rx="4" fill="{col}" opacity="0.12"/>')
    rows2 += (f'<text x="{COL_DROSO + COL_W/2:.0f}" y="{y_center - 5:.1f}" '
              f'text-anchor="middle" font-size="12" font-weight="700" fill="{col}" '
              f'font-style="italic">{g["drosophila_ortholog"]}</text>')
    rows2 += (f'<text x="{COL_DROSO + COL_W/2:.0f}" y="{y_center + 11:.1f}" '
              f'text-anchor="middle" font-size="8" fill="#666">'
              f'{g["journal"]} {g["year"]}</text>')

    # Arrow 2
    arrow2_x1 = COL_DROSO + COL_W + 4
    arrow2_x2 = COL_PHENO - 4
    rows2 += (f'<line x1="{arrow2_x1}" y1="{y_center:.1f}" x2="{arrow2_x2}" y2="{y_center:.1f}" '
              f'stroke="{col}" stroke-width="1.5"/>')
    rows2 += (f'<text x="{(arrow2_x1+arrow2_x2)/2:.0f}" y="{y_center - 5:.1f}" '
              f'text-anchor="middle" font-size="7.5" fill="{col}">RNAi / mutant</text>')

    # Phenotype box
    rows2 += (f'<rect x="{COL_PHENO}" y="{y_top:.1f}" width="{COL_W}" '
              f'height="{row_h-12:.1f}" rx="4" fill="{col}" opacity="0.12"/>')
    pheno = g["fly_phenotype"]
    rows2 += (f'<text x="{COL_PHENO + COL_W/2:.0f}" y="{y_center + 4:.1f}" '
              f'text-anchor="middle" font-size="9" fill="{col}" font-weight="600">'
              f'{pheno[:24]}</text>')

    # Mechanism (truncated)
    mech = g["mechanism"][:35]
    rows2 += (f'<text x="{COL_MECH + 8}" y="{y_center + 4:.1f}" '
              f'font-size="8.5" fill="#444">{mech}</text>')

# Defs for arrowheads (not used in SVG without marker support, but add lines)
defs = ""
svg2 = f"""<svg viewBox="0 0 {FW2} {FH2}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW2//2}" y="18" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Kayser Lab Cross-Species Validation Pipeline
  </text>
  <text x="{FW2//2}" y="34" text-anchor="middle" font-size="10" fill="#666">
    Human GWAS candidate gene → Drosophila ortholog knockdown → sleep phenotype → circuit mechanism
  </text>
  {rows2}
  <text x="{FW2//2}" y="{FH2 - 8}" text-anchor="middle" font-size="8.5" fill="#aaa">
    Data: Mace et al. 2026 Nat Commun · Mace et al. 2026 biorXiv · Rodriguez et al. 2026 biorXiv · Coll-Tané et al. 2026 J Clin Invest
  </text>
</svg>"""

with open(os.path.join(OUT, "validation_pipeline.svg"), "w") as f:
    f.write(svg2)
print("Wrote validation_pipeline.svg")
