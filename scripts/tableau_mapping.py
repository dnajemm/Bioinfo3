import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
"""Création d'un tableau récapitulatif des statistiques de mapping (taux de mapping et profondeur moyenne) pour toutes les souches.
Le tableau est sauvegardé en PNG dans le dossier results_mapping.""""

# Chemins 
INPUT = "/data/projet3/Output/results_mapping/recap_mapping.csv"
OUTDIR = Path("/data/projet3/Output/results_mapping")

# Charger le CSV
df = pd.read_csv(INPUT)

# Calculer les moyennes sur toutes les souches
mean_mapping = df["Taux_mapping(%)"].mean()
mean_depth = df["Profondeur_moyenne"].mean()

# Créer un petit dataframe récapitulatif des moyennes
df_mean = pd.DataFrame({
    "Paramètre": ["Taux de mapping moyen (%)", "Profondeur moyenne (X)"],
    "Valeur moyenne": [round(mean_mapping, 2), round(mean_depth, 2)]})

# Créer la figure
fig, ax = plt.subplots(figsize=(6, 2))
ax.axis("off")

# Créer le tableau avec uniquement les moyennes
table = ax.table(
    cellText=df_mean.values,
    colLabels=df_mean.columns,
    loc="center")

table.auto_set_font_size(False)
table.set_fontsize(10)
table.auto_set_column_width(col=list(range(len(df_mean.columns))))

# Couleur du header 
for (row, col), cell in table.get_celld().items():
    if row == 0:  # ligne d'en-tête
        cell.set_facecolor("#bcd4e6")  # bleu clair
        cell.set_text_props(weight='bold', color='black')

# Sauvegarder dans le dossier OUTDIR 
outfile = OUTDIR / "mapping_summary.png"
plt.tight_layout()
plt.savefig(outfile, dpi=300, bbox_inches="tight")
