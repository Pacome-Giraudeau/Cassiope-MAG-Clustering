import os
import sys
import numpy as np

def analyse_distances():
    # Récupération du chemin du dossier contenant le script
    folder = os.path.dirname(os.path.abspath(sys.argv[0]))
    filepath = os.path.join(folder, "data", "first_set", "distance_data_analyse.txt")

    # Lecture des distances depuis le fichier
    distances = []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                try:
                    d = float(parts[2])
                    distances.append(d)
                except ValueError:
                    continue  # Ignore les lignes mal formées

    if not distances:
        print("Aucune distance valide trouvée.")
        return

    # Calcul de la moyenne
    moyenne = np.mean(distances)

    # Calcul du 80e percentile (seuil pour lequel 80% des distances sont ≤ à ce seuil)
    seuil_80 = np.percentile(distances, 80)

    # Résultats
    print(f"Nombre total de distances : {len(distances)}")
    print(f"Distance moyenne : {moyenne:.4f}")
    print(f"Seuil 80% : {seuil_80:.4f} (80% des distances ≤ à ce seuil)")

if __name__ == "__main__":
    analyse_distances()
