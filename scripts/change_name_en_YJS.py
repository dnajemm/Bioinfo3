import os 
import pandas as pd
#Ce code permet de renommer les fichiers fastq en leur nom de souche d'après le dictionnaire name_dict contenant les anciens et nouveaux noms selon le fichier excel.
#On commence par créer le dictionnaire de correspondance entre les anciens et nouveaux noms de fichiers fastq.
#Il lit les noms depuis un fichier Excel et génère les paires de noms avec les suffixes _1 et _2.

# Charger le fichier Excel
df = pd.read_excel("/data/projet3/Input/tableauSouches.xlsx")

# Crée le dictionnaire de base avec les colonnes 0 et 7 du DataFrame qui contiennent les anciens et nouveaux noms.
base_dict = dict(zip(df.iloc[:, 0], df.iloc[:, 7]))
# Nouveau dictionnaire avec les suffixes _1 et _2
name_dict = {}
for yjs, err in base_dict.items(): # yjs = nouveau nom, err = ancien nom
    name_dict[f"{err}_1.fastq"] = f"{yjs}_1.fastq"
    name_dict[f"{err}_2.fastq"] = f"{yjs}_2.fastq"

# Dossier contenant les fichiers à renommer
folder = "/data/projet3/Input/reads_paired_end/reads/"

for old_name, new_name in name_dict.items(): 
    old_path = os.path.join(folder, old_name) # chemin complet de l'ancien fichier
    new_path = os.path.join(folder, new_name) # chemin complet du nouveau fichier
     # renommer le fichier s'il existe
    if os.path.exists(old_path):
        os.rename(old_path, new_path)
        print(f" {old_name} → {new_name}")
    else: # le fichier n'existe pas
        print(f" {old_name} introuvable")
