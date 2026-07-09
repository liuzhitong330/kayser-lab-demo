"""
fetch_data.py — Kayser Lab Demo (Matthew Kayser, UCSF)
Downloads sleep trait GWAS landscape from the GWAS Catalog (EBI) and combines with
manually curated data from the Kayser lab's published cross-species validation studies.

The Kayser lab pipeline: human sleep GWAS → Drosophila functional screen → circuit mechanism.
"""

import urllib.request, csv, io, os
from collections import defaultdict

OUT = os.path.dirname(__file__)

# ── 1. GWAS Catalog sleep trait landscape ────────────────────────────────────
print("Downloading GWAS Catalog studies …")
url = "https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative"
with urllib.request.urlopen(url, timeout=60) as r:
    content = r.read().decode("utf-8")
print(f"  Downloaded {len(content):,} bytes")

SLEEP_KW = [
    "insomnia", "chronotype", "narcolep", "excessive daytime sleepiness",
    "sleep duration", "sleep latency", "sleep quality", "hypersomnia",
    "morning person", "evening person", "ease of getting up"
]
TRAIT_CANON = {
    "insomnia measurement": "insomnia",
    "insomnia": "insomnia",
    "chronotype measurement": "chronotype",
    "sleep duration": "sleep duration",
    "narcolepsy-cataplexy syndrome": "narcolepsy",
    "narcolepsy without cataplexy": "narcolepsy",
    "narcolepsy": "narcolepsy",
    "excessive daytime sleepiness measurement": "excessive daytime sleepiness",
    "hypersomnia": "hypersomnia",
    "sleep latency": "sleep latency",
    "sleep quality": "sleep quality",
    "ease of getting up in the morning": "ease of getting up",
}

reader = csv.DictReader(io.StringIO(content), delimiter="\t")
trait_data = defaultdict(lambda: {"studies": 0, "assocs": 0})
for row in reader:
    trait_raw = row.get("MAPPED_TRAIT", "").lower().split(",")[0].strip()
    if not any(kw in trait_raw for kw in SLEEP_KW):
        continue
    canon = TRAIT_CANON.get(trait_raw, trait_raw[:40])
    trait_data[canon]["studies"] += 1
    try:
        trait_data[canon]["assocs"] += int(row.get("ASSOCIATION COUNT", 0) or 0)
    except ValueError:
        pass

# Keep main traits only
KEEP = ["insomnia", "chronotype", "sleep duration", "narcolepsy",
        "excessive daytime sleepiness", "hypersomnia", "sleep latency", "sleep quality"]
trait_rows = [(t, trait_data[t]["assocs"], trait_data[t]["studies"])
              for t in KEEP if t in trait_data]
trait_rows.sort(key=lambda x: -x[1])

with open(os.path.join(OUT, "sleep_gwas_traits.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["trait", "gwas_loci", "n_studies"])
    for t, a, s in trait_rows:
        w.writerow([t, a, s])
print("Wrote sleep_gwas_traits.tsv")
print()
for t, a, s in trait_rows:
    print(f"  {t:40s}: {a:5d} loci, {s} studies")

# ── 2. Kayser lab validated genes ─────────────────────────────────────────────
# Manually curated from published papers (verified July 2026):
# - Mace et al. 2026 Nat Commun: CADM2/beat-Ia, excessive daytime sleepiness
# - Mace et al. 2026 biorXiv: EGR2/stripe, idiopathic hypersomnia (cross-trait GWAS)
# - Rodriguez et al. 2026 biorXiv: SIK3/Sik3, sleep disruption in TDP-43 proteinopathy
# - Coll-Tané et al. 2026 J Clin Invest: FOXP1/FoxP, insomnia / FOXP1 syndrome
VALIDATED_GENES = [
    # human_gene, droso_ortholog, trait, phenotype_dir, mechanism, paper_year, journal
    ("CADM2",  "beat-Ia", "excessive daytime sleepiness", "increased sleep",
     "NPF neurite defect → arousal circuit failure",    2026, "Nat Commun"),
    ("EGR2",   "stripe",  "hypersomnia (sub-threshold)", "increased sleep + consolidation",
     "EGR2 as distal effector at ADO-EGR2 locus",      2026, "biorXiv"),
    ("SIK3",   "Sik3",    "sleep disruption (TDP-43)",   "loss rescues sleep",
     "Peripheral metabolic dysfunction via SIK3",       2026, "biorXiv"),
    ("FOXP1",  "FoxP",    "insomnia",                    "fragmented/reduced sleep",
     "Developmental peptidergic signaling (PDF neurons)",2026, "J Clin Invest"),
]

with open(os.path.join(OUT, "validated_genes.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["human_gene", "drosophila_ortholog", "sleep_trait", "fly_phenotype",
                "mechanism", "year", "journal"])
    for row in VALIDATED_GENES:
        w.writerow(row)
print("\nWrote validated_genes.tsv")
