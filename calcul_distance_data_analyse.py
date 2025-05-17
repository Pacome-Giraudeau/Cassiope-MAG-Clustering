import os
import sys
import numpy as np

def analyse_distances():
    # Récupération du chemin du dossier contenant le script
    folder = os.path.dirname(os.path.abspath(sys.argv[0]))
    filepath = os.path.join(folder, "data", "first_set", "distance_data_analyse.txt")
    seuil_file = os.path.join(folder, "data", "first_set", "seuil.txt")

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

    # Calcul du 80e percentile
    seuil_80 = np.percentile(distances, 80)

    # Résumé
    print(f"Nombre total de distances : {len(distances)}")
    print(f"Distance moyenne : {moyenne:.4f}")
    print(f"Seuil 80% : {seuil_80:.4f} (80% des distances ≤ à ce seuil)")

    # Écriture des percentiles 95, 90, ..., 10 dans le fichier seuil.txt
    with open(seuil_file, "w") as f:
        for p in range(95, 9, -5):  # de 95 à 10 inclus, par pas de -5
            val = np.percentile(distances, p)
            f.write(f"{p}\t{val:.4f}\n")

    print(f"\nFichier seuil.txt écrit avec les seuils de 95% à 10% (par pas de 5%).")

    


if __name__ == "__main__":
    analyse_distances()
