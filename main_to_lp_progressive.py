import clingo
import sys
import os
import numpy as np
from scipy.spatial.distance import euclidean
import subprocess
import random
import glob

###### PARAMÈTRES

nb_clusters = 2
step_size = 100  # nombre de contigs par itération
max_total_contigs = 300  # nombre total de contigs à traiter
classified_contigs = {}  # contig -> cluster

folder = os.path.dirname(os.path.abspath(sys.argv[0]))
folder_data = "data"
folder_data_set = "first_set"
folder_data_working = os.path.join(folder, folder_data, folder_data_set)

def copy_random_lines_with_scg(folder_data_working, nb_contigs, scg_file="contig_scg_details.tsv", output_suffix="_subset"):
    # Chemins des fichiers source
    file1 = os.path.join(folder_data_working, f"kmer_contigs.clean.txt")
    file2 = os.path.join(folder_data_working, f"mean_coverage_Q2Q3_contigs.txt")
    scg_path = os.path.join(folder_data_working, scg_file)

    # Fichiers de sortie
    output1 = os.path.join(folder_data_working, f"{nb_contigs}_kmer_contigs{output_suffix}.txt")
    output2 = os.path.join(folder_data_working, f"{nb_contigs}_mean_coverage_Q2Q3{output_suffix}.txt")

    # Charger les contigs SCG (premier mot de chaque ligne)
    with open(scg_path, 'r') as f_scg:
        scg_contigs = {line.strip().split()[0] for line in f_scg if line.strip()}

    print(f"Contigs SCG trouvés : {len(scg_contigs)}")

    # Lire les deux fichiers et créer des dictionnaires {contig_name: line}
    def read_file_to_dict(file_path, separator):
        with open(file_path, 'r') as f:
            lines = f.readlines()[1:]
        return {line.split(separator)[0]: line for line in lines}

        
    def read_file_to_dict_coverage(file_path, separator):
        with open(file_path, 'r') as f:
            lines = f.readlines()
        return {(line.split(separator)[0], line.split(separator)[1]): line for line in lines}

    contig_to_line1 = read_file_to_dict(file1, "	")
    contig_to_line2 = read_file_to_dict_coverage(file2, "	") 

    # Vérifier que les contigs sont cohérents entre les deux fichiers
    contigs_common = set(contig_to_line1.keys())

    # Extraire les lignes SCG
    scg_contigs_present = scg_contigs & contigs_common
    scg_lines1 = [contig_to_line1[contig] for contig in scg_contigs_present]

    scg_lines2 = [contig_to_line2[contig+"_split_00001", echantillon] for contig in scg_contigs_present for echantillon in [
        "S_125SUR1QQSS11",
        "S_133DCM1SSUU11",
        "S_135SUR1MMQQ11",
        "S_145SUR1SSUU11", 
        "S_152SUR1SSUU11",
        "S_66SUR1SSUU11",
        "S_81SUR01SSUU11"
    ]]

    # Contigs restants (non-SCG)
    remaining_contigs = list(contigs_common - scg_contigs_present)

    if nb_contigs < len(scg_contigs_present):
        raise ValueError(f"Le nombre de contigs SCG ({len(scg_contigs_present)}) est supérieur à {nb_contigs}.")

    if nb_contigs - len(scg_contigs_present) > len(remaining_contigs):
        raise ValueError(f"Pas assez de contigs hors SCG pour sélectionner {nb_contigs - len(scg_contigs_present)} contigs.")

    # Sélection aléatoire de contigs non-SCG
    random_selected_contigs = random.sample(remaining_contigs, nb_contigs - len(scg_contigs_present))

    # Construire les lignes finales
    final_lines1 = scg_lines1 + [contig_to_line1[contig] for contig in random_selected_contigs]
    final_lines2 = scg_lines2 + [contig_to_line2[contig+"_split_00001", echantillon] for contig in random_selected_contigs for echantillon in [
        "S_125SUR1QQSS11",
        "S_133DCM1SSUU11",
        "S_135SUR1MMQQ11",
        "S_145SUR1SSUU11", 
        "S_152SUR1SSUU11",
        "S_66SUR1SSUU11",
        "S_81SUR01SSUU11"
    ]]

    # Écrire les fichiers de sortie
    with open(output1, 'w') as f_out1:
        f_out1.writelines(final_lines1)

    with open(output2, 'w') as f_out2:
        f_out2.writelines(final_lines2)

    print(f"Fichiers écrits :\n- {output1}\n- {output2}")

copy_random_lines_with_scg(folder_data_working=folder_data_working, nb_contigs=step_size)

""" Description du chemin vers scg.expected qui contient 
la liste des gènes marqueurs bactériens attendus
dans un génome"""
tokens = dict() # Association gène marqueur, groupe taxonomique
tokens_files = [
    os.path.join(folder, folder_data, folder_data_set, "scg.expected")
]

""" Description du chemin vers mean_coverage_Q2Q3_contigs qui contient 
les informations d'abondance par échantillon de chacun des contigs """
contigs_coverage = dict()
contigs_coverage_files = os.path.join(folder, folder_data, folder_data_set, f"{step_size}_mean_coverage_Q2Q3_subset.txt")

""" Description du chemin vers contig_scg_details qui contient 
les gènes marqueurs identifiés pour chaque contig et leurs affiliations taxonomiques """
contigs_tokens = dict() # Association contigs, gène marqeurs, groupe taxonomique détaillé
contigs_tokens_files = [
    os.path.join(folder, folder_data, folder_data_set, "contig_scg_details.tsv")
]

""" Description du chemin vers kmer_contigs.clean_very_short qui contient 
une version raccourcie de la table des fréquences en
tétranucléotides de chacun des contigs """
contigs_kmere = dict() # Décompte des bases tétranucléotide de chaque contig
contigs_kmere_files = [
    os.path.join(folder, folder_data, folder_data_set, f"{step_size}_kmer_contigs_subset.txt")
]

contigs = list() # Initialisation de la liste des contigs



####### LECTURE DES DONN2ES

def read_tokens():
    """ Lecture du fichier scg.expected et extraction des gènes marqueurs bactériens.
    Input : Fichier scg.expected
    Output : Dictionnaire {gène marqueur : groupe taxonomique}"""
    for t_file in tokens_files:
        with open(t_file, "r", encoding="utf-8") as file:
            for line in file:
                try:
                    tokens[line.split()[0]].append(line.split()[1]) # ajout des gènes marqueurs extrait du fichier scg.expected
                except KeyError:
                    tokens[line.split()[0]] = [line.split()[1]]
    return tokens
                
'''def read_contig_scg():
    print("Début de read_contig_scg")
    for scg_file in contigs_tokens_files:
        with open(scg_file, "r", encoding="utf-8") as file:
            for line in file:
                data = line.split()
                contig_name = data[0]  
                gene = data[1]
                taxonomie = data[4:11]
                
                # On cherche si le contig est déjà dans la matrice
                contig_found = False
                for entry in contigs_tokens:
                    if entry[0] == contig_name: # si oui, on ajoute le gène à la liste des gènes marqueurs du contigs
                        entry[1].append(gene)
                        contig_found = True
                        break

                # Si le contig n'est pas dans la matrice, on l'ajoute
                if not contig_found:
                    contigs_tokens.append([contig_name, [gene], taxonomie]) # si non, on ajoute l'association contig, gène marqueur, taxinomie détaillée

    # Affiche chaque contig avec ses gènes et sa taxonomie pour vérification
    for entry in contigs_tokens:
        print(f"Contig: {entry[0]} | Gènes: {entry[1]} | Taxonomie: {entry[2]}") # On affiche le résultat
    print(f"Fin de read_contig_scg(), contigs_tokens contient {len(contigs_tokens)} entrées")'''

#SARRA : nouveau read_contig_scg()
def read_contig_scg():
    """Lecture du fichier contig_scg_details.tsv et extraction des données taxonomiques + gènes marqueurs.
    Input : Fichier contig_scg_details.tsv
    Output : Dictionnaire {contig : ([gènes marqueurs], [taxonomie])}"""
    global contigs_tokens  # Dictionnaire contenant les contigs avec gènes marqueurs et taxonomie
    
    for scg_file in contigs_tokens_files:
        print(f"Lecture du fichier : {scg_file}")  
        with open(scg_file, "r", encoding="utf-8") as file:
            next(file)  # Sauter l'en-tête
            for line in file:
                data = line.strip().split("\t")  # Séparer par tabulation
                contig_name = data[0]  # Colonne "contig"
                gene = data[1]  # Colonne "gene_name"
                taxonomie = data[5:11]  # Colonnes t_domain → t_genus

                # Vérifier si le contig est déjà enregistré
                if contig_name in contigs_tokens:
                    contigs_tokens[contig_name][0].append(gene)  # Ajouter le gène à la liste existante
                else:
                    contigs_tokens[contig_name] = ([gene], taxonomie)  # Créer une nouvelle entrée (liste de gènes, taxonomie)

    # Vérification : Affichage des 5 premiers contigs stockés
    print("\n Vérification des données stockées :")
    for i, (contig, (genes, taxo)) in enumerate(contigs_tokens.items()):
        # print(f"{i+1}. Contig: {contig}")
        # print(f"   Gènes marqueurs: {genes}")
        # print(f"   Taxonomie: {taxo}\n")
        if i == 4:  # Limite à 5 affichages
            break
    return contigs_tokens

def read_contig_kmere():
    """Lecture du fichier kmer_contigs.clean.txt et extraction des données tétranucléotides.
    Input : Fichier kmer_contigs.clean
    Output : Dictionnaire {contig : [fréquences tétranucléotides]}, Liste [contigs]"""
    liste_contig = []  # Liste pour stocker les contigs dans l'ordre de lecture
    for kmere_file in contigs_kmere_files:
        with open(kmere_file, "r", encoding="utf-8") as file:
            for line in file:
                if line.split()[0] != "contig":
                    contig_name = line.split()[0]
                    liste_contig.append(contig_name)  # Ajout du contig à la liste
                    contigs_kmere[contig_name] = np.array(line.split()[1:]).astype(int)  # Stockage des kmères sous forme d'entiers
    # print(liste_contig[:5])  # Affiche les 5 premiers contigs
    return contigs_kmere, liste_contig

def read_contig_coverage():
    """Lecture du fichier mean_coverage_Q2Q3_contigs_2.txt et extraction des données d'abondance.
    Input : Fichier mean_coverage_Q2Q3_contigs_2.txt
    Output : Matrice (contigs x stations), Liste [contigs]"""
    stations = []  # Utilisation d'une liste pour collecter toutes les stations
    data_dict = {}  # Dictionnaire temporaire pour stocker les valeurs
    liste_contig = []

    # Lecture du fichier pour récupérer les contigs et les stations
    with open(contigs_coverage_files , "r", encoding="utf-8") as file:
        for line in file:
            parts = line.strip().split("\t")
            contig_split, station, coverage = parts
            contig = contig_split.split("_split_")[0]
            coverage = float(coverage.replace(",", "."))  # Convertit en float

            if contig not in data_dict:
                data_dict[contig] = {}
                liste_contig.append(contig)  # Ajoute le contig dans l'ordre de lecture
            data_dict[contig][station] = coverage
            if station not in stations:
                stations.append(station)

    # Construction de la matrice de taille (nb_contigs x nb_stations)
    contigs_coverage = np.zeros((len(liste_contig), len(stations)))

    # Remplissage de la matrice
    index = 0
    for contig in liste_contig:
        for station in data_dict[contig]:
            j = stations.index(station)
            contigs_coverage[index, j] = data_dict[contig][station]
        index += 1

    # Affichage d'un extrait de la matrice
    # print("Matrice des contigs (shape):", contigs_coverage.shape)
    # print(contigs_coverage[:5, :5])  # Affiche les 5 premières lignes et colonnes
    return contigs_coverage, liste_contig


def calculate_distance_contigs_naif():
    """Calcul de la distance euclidienne entre les contigs à partir des données d'abondance et de tétranucléotides.
    Création d'un vecteur combiné pour chaque contig et calcul de la distance euclidienne entre chaque paire de contigs.
    Input : contigs_coverage, contigs_kmere
    Output : Matrice des distances, Dictionnaire {contig : index}"""
    # Chargement des données
    contigs_coverage, liste_coverage = read_contig_coverage()
    contigs_kmere, liste_kmere = read_contig_kmere()

    # Identification des contigs en commun
    contigs = list(set(liste_coverage) & set(liste_kmere))
    n = nb_contigs
    
    # Initialisation de la matrice de distances
    contigs_dist = np.zeros((n, n))

    # Création d'un dictionnaire pour indexer les contigs
    contig_to_index = {contig: idx for idx, contig in enumerate(contigs)}

    for i in range(n):
        c1 = contigs[i]
        # print(i + 1, "/", n)

        # Construction du vecteur combiné pour c1
        vecteur_c1 = np.concatenate((contigs_coverage[liste_coverage.index(c1)], contigs_kmere[c1]))

        for j in range(i, n):  # Évite les doublons (distance symétrique)
            c2 = contigs[j]
            vecteur_c2 = np.concatenate((contigs_coverage[liste_coverage.index(c2)], contigs_kmere[c2]))

            # Calcul de la distance euclidienne
            distance = euclidean(vecteur_c1, vecteur_c2)
            contigs_dist[i, j] = distance
            contigs_dist[j, i] = distance  # Symétrie

    print("Matrice des distances (extrait 5x5):")
    print(contigs_dist[:5, :5])

    return contigs_dist, contigs


def calculate_distance_contigs():

    """Calcul de la distance euclidienne entre les contigs à partir des données d'abondance et de tétranucléotides.
    Calcul de la distance euclidienne entre les vecteurs d'abondance puis de tétranucléotides de chaque paire de contigs
    puis moyenne avec pondération.
    Input : contigs_coverage, contigs_kmere
    Output : Matrice des distances, Dictionnaire {contig : index}"""

    ponderation = 300
    # Chargement des données
    contigs_coverage, liste_coverage = read_contig_coverage()
    contigs_kmere, liste_kmere = read_contig_kmere()

    # Identification des contigs en commun
    contigs_en_commun = list(set(liste_coverage) & set(liste_kmere))
    # print("Liste des contigs en commun :", contigs_en_commun)
    for contig in contigs_en_commun :
        contigs.append(contig)
    n = nb_contigs
    
    # Initialisation de la matrice de distances
    contigs_dist = np.zeros((n, n))

    # Création d'un dictionnaire pour indexer les contigs
    contig_to_index = {contig: idx for idx, contig in enumerate(contigs)}

    # Préparation du fichier des distances
    output_path = os.path.join(folder_data_working, "distance_data_analyse.txt")


    # Calcul de la matrice des distances
    with open(output_path, "w") as f:
        # Calcul de la matrice des distances et écriture dans le fichier
        for i in range(n):
            c1 = contigs[i]
            print(f"{i + 1} / {n}")
            coverage_c1 = contigs_coverage[liste_coverage.index(c1)]
            kmere_c1 = contigs_kmere[c1]

            for j in range(i, n):  # Évite les doublons (distance symétrique)
                c2 = contigs[j]
                coverage_c2 = contigs_coverage[liste_coverage.index(c2)]
                kmere_c2 = contigs_kmere[c2]

                # Calcul de la distance euclidienne pondérée
                distance_coverage = euclidean(coverage_c1, coverage_c2)
                distance_kmere = euclidean(kmere_c1, kmere_c2)
                distance = (distance_coverage * ponderation + distance_kmere) / (ponderation + 1)
                contigs_dist[i, j] = distance

                # Écriture dans le fichier si i != j
                if i != j:
                    f.write(f"{c1}\t{c2}\t{distance:.6f}\n")
            
    # print("Matrice des distances (extrait 5x5):")
    # print(contigs_dist[:5, :5])

    return contigs_dist, contigs



def get_distance_contigs(c1, c2, contigs_dist):
    """
    c1 et c2 deux entiers
    Renvoie la distance entre c1 et c2
    """
    return int(contigs_dist[min(c1, c2), max(c1, c2)])

def get_nb_contigs():
    """
    renvoie le nombre de contigs considérés
    """
    nb = 0
    for kmere_file in contigs_kmere_files:
        with open(kmere_file, "r", encoding="utf-8") as file:
            for line in file:
               nb+=1 
    return nb 
def get_contigs():
    """
    renvoie la liste des contigs
    """
    liste_contig = []  # Liste pour stocker les contigs dans l'ordre de lecture
    for kmere_file in contigs_kmere_files:
        with open(kmere_file, "r", encoding="utf-8") as file:
            for line in file:
                if line.split()[0] != "contig":
                    contig_name = line.split()[0]
                    liste_contig.append(contig_name)  # Ajout du contig à la liste
    return liste_contig

    ##################  
    ##################  
######### CLINGO ######### 
    ##################  
    ##################  


def add_clingo_contigs():
    """Ajoute les contigs à Clingo.
    Input : Liste des contigs
    Output : Règles Clingo contig("contig1")."""

    print("Cration des contigs...")
    
    file_path = os.path.join(folder_data_working, "programme_lp", str(nb_contigs) + "_contigs.lp")
    if not os.path.exists(file_path):
        contigs = get_contigs()

        print("Vérification de la liste des contigs :", contigs)
        with open(file_path, "w") as f:
            for c1 in contigs:
                print(f"Ajout du contig {c1} dans le fichier contig.lp")
                f.write(f'contig("{c1}").\n')
    # ctl.load(file_path)
    print("Création des contigs effectué !\n")
    return file_path


        

def add_clingo_restriction_tokens():
    """Ajoute les distances entre les contigs à Clingo.
    Input : Matrice des distances, Dictionnaire {contig : index}
    Output : Règles Clingo distance("contig1", "contig2", 5)."""
    
    print("Création des restrictions via tokens...")

    
    file_path = os.path.join(folder_data_working, "programme_lp", str(nb_contigs) + "_restrictions_via_tokens.lp")
    if not os.path.exists(file_path):
        contigs_dist, contigs = calculate_distance_contigs()
        
        with open(file_path, "w") as f:
            for c1 in range(nb_contigs):
                print(f"[Ajout des distances restrictions [{c1} / {nb_contigs-1}]]")

                # Calcule toutes les distances de c1 à chaque c2
                distances = []
                for c2 in contigs_tokens.keys():
                    distance = get_distance_contigs(c1, contigs.index(c2), contigs_dist)
                    distances.append((distance, c2))
                
                # Trie les c2 par distance croissante
                distances.sort()
                # Prend la moitié des c2 les plus proches
                top_half = distances[:len(distances) // 2]
                
                # Écrit la contrainte seulement pour cette moitié
                for _, c2 in top_half:
                    f.write(f':- assigned("{contigs[c1]}", C), assigned("{c2}", C), cluster(C).\n')

    print("Création des restrictions effectuées !\n")
    
    # ctl.load(file_path)
    return file_path

        

def add_clingo_distance():
    """Ajoute les distances entre les contigs à Clingo.
    Input : Matrice des distances, Dictionnaire {contig : index}
    Output : Règles Clingo distance("contig1", "contig2", 5)."""
    
    print("Création des distances...")
    
    file_path = os.path.join(folder_data_working, "programme_lp", str(nb_contigs) + "_distances.lp")
    if not os.path.exists(file_path):
        contigs_dist, contigs = calculate_distance_contigs()
        
        with open(file_path, "w") as f:
            for c1 in range(0, nb_contigs):
                
                print("[Ajout des distances entre contigs [", c1, "/", nb_contigs,"]")
                for c2 in range(nb_contigs):
                    f.write(f'distance("{contigs[c1]}","{contigs[c2]}", {get_distance_contigs(c1, c2, contigs_dist)}).\n')
    
    print("Création des distances effectué !\n")
    
    # ctl.load(file_path)
    return file_path


def add_clingo_scg_and_taxonomy():
    """Ajoute les gènes marqueurs et la taxonomie des contigs à Clingo.
    Input : Dictionnaire {contig : ([gènes marqueurs], [taxonomie])}
    Output : Règles Clingo scg("contig1", "gene1")."""

    print("Création des contigs marqueurs et de leurs taxonomie...")
    
    file_path = os.path.join(folder_data_working, "programme_lp", str(nb_contigs) + "_scg_&_taxonomy.lp")
    if not os.path.exists(file_path):
        with open(file_path, "w") as f:
            
            for contig, (genes, taxo) in contigs_tokens.items():
                for gene in genes:
                    #print(f'Ajout Clingo : scg("{contig}", "{gene}")')
                    f.write(f'scg("{contig}", "{gene}").\n')

                for rank, taxon in enumerate(taxo, start=1):
                    #print(f'Ajout Clingo : taxo("{contig}", {rank}, "{taxon}")')
                    f.write(f'taxo("{contig}", {rank}, "{taxon}").\n')
    print("Création des taxonomies effectué !\n")
    
    # ctl.load(file_path)
    
    return file_path

def add_clingo_clusters():
    """Ajoute les clusters à Clingo.
    Input : Nombre de clusters  
    Output : Règles Clingo cluster(1)."""


    file_path = os.path.join(folder_data_working, "programme_lp", "clusters.lp")

    print("Création des clusters...")
    with open(file_path, "w") as f:
        for cluster in range(nb_clusters):
            f.write(f'cluster({cluster}).\n') 

    print("Création des clusters effectué !\n")
    # ctl.load(file_path)
    return file_path

def add_clingo_rules():
    """Ajoute les clusters à Clingo.
    Input : Nombre de clusters  
    Output : Règles Clingo cluster(1)."""


    file_path = os.path.join(folder, "rules.lp")

    # ctl.load(file_path)
    return file_path

#LISTE DES FICHIERS RULES
#rules_files = glob.glob(os.path.join(folder, "rules", "**", "rules_*.lp"), recursive=True)
#rules_files = glob.glob(os.path.join(folder, "rules", "test", "rules_*.lp"))
#rules_files = ["rules.lp"]
rules_files = glob.glob(os.path.join(folder, "rules", "progressive", "**", "rules_*.lp"), recursive=True)
print(f"Fichiers rules détectés :")
for rule in rules_files:
    print(" -", rule)

# FICHIER FIXED 
def write_fixed_assignments(fixed_path):
    with open(fixed_path, "w") as f:
        if classified_contigs:
            for contig, cluster in classified_contigs.items():
                f.write(f'fixed("{contig}", {cluster}).\n')
            f.write("\nassigned(C,B) :- fixed(C,B).\n")
            f.write(":- assigned(C,B1), fixed(C,B2), B1 != B2.\n")


# ON_MODEL 
def on_model(model, rules_name):
    global model_count
    model_count += 1

    atoms = model.symbols(shown=True)
    assigned_atoms = [atom for atom in atoms if atom.name == "assigned" and len(atom.arguments) == 2]

    # Mise à jour de classified_contigs
    for atom in assigned_atoms:
        contig = str(atom.arguments[0]).strip('"')
        cluster = int(str(atom.arguments[1]))
        classified_contigs[contig] = cluster

    unique_clusters = set(cluster for _, cluster in classified_contigs.items())
    count_clusters = len(unique_clusters)
    count_assigned = len(classified_contigs)

    # Affichage des résultats
    print(f"\nRésultat pour {rules_name} :")
    print(f"Nombre de clusters distincts : {count_clusters}")
    print(f"Nombre total de contigs assignés (cumulés) : {count_assigned}")
    print("Liste des assignations dans ce modèle :")
    for atom in assigned_atoms:
        print(str(atom))

    print(f"\nOptimalité prouvée ? {model.optimality_proven}")
    print("Fin de l'affichage des clusters\n")

    # Sauvegarde cumulative
    clusters_output_file = os.path.join(folder_data_working, f"{len(classified_contigs)}_clusters_{rules_name}.txt")
    with open(clusters_output_file, "w") as f:
        for contig, cluster in classified_contigs.items():
            f.write(f"{contig}\tbin_{cluster}\n")

    if model_count >= 5:
        return False


# INIT_CLINGO
def init_clingo(rules_file_path, rules_name):
    file_clusters = add_clingo_clusters()
    file_contigs = add_clingo_contigs()
    file_distances = add_clingo_distance()
    file_scg_taxo = add_clingo_scg_and_taxonomy()
    fixed_file = os.path.join(folder_data_working, "programme_lp", "fixed.lp")

    file_output = os.path.join(folder, f"output_{rules_name}.lp")
    
    print("Merging les fichiers suivants dans :", file_output)
    print(" -", file_clusters)
    print(" -", file_contigs)
    print(" -", file_distances)
    print(" -", file_scg_taxo)
    print(" -", fixed_file)
    print(" -", rules_file_path)

    with open(file_output, "w") as f_out:
        with open(file_clusters, "r") as f1:
            f_out.write(f1.read() + "\n")
        with open(file_contigs, "r") as f2:
            f_out.write(f2.read() + "\n")
        with open(file_distances, "r") as f3:
            f_out.write(f3.read() + "\n")
        with open(file_scg_taxo, "r") as f4:
            f_out.write(f4.read() + "\n")
        with open(fixed_file, "r") as f5:
            f_out.write(f5.read() + "\n")
        with open(rules_file_path, "r") as f6:
            f_out.write(f6.read() + "\n")

    print("Création du fichier programme OK pour :", rules_name)
    return file_output

def load_all_contigs_source():
    source_file = os.path.join(folder_data_working, "kmer_contigs.clean.txt")
    with open(source_file, "r") as f:
        return [line.split("\t")[0] for line in f.readlines()[1:]]

all_contigs_source = load_all_contigs_source()

# MAIN
def __main__():
    global model_count, nb_contigs, contigs_kmere_files, contigs_coverage_files

    step_id = 1
    nb_total = 0

    while nb_total < max_total_contigs:
        print(f"\n======== ÉTAPE {step_id} ========")
        
        available_contigs = [c for c in all_contigs_source if c not in classified_contigs]
        print(f"Contigs disponibles restants : {len(available_contigs)}")

        if not available_contigs:
            print("Plus de contigs disponibles.")
            break

        step_contigs = random.sample(available_contigs, min(step_size, max_total_contigs - nb_total))
        nb_contigs = len(step_contigs)

        copy_random_lines_with_scg(folder_data_working, nb_contigs)

        contigs_kmere_files = [os.path.join(folder_data_working, f"{nb_contigs}_kmer_contigs_subset.txt")]
        contigs_coverage_files = os.path.join(folder_data_working, f"{nb_contigs}_mean_coverage_Q2Q3_subset.txt")

        read_contig_kmere()
        read_contig_scg()

        fixed_path = os.path.join(folder_data_working, "programme_lp", "fixed.lp")
        write_fixed_assignments(fixed_path)

        for rules_path in rules_files:
            model_count = 0
            rules_name = f"{os.path.splitext(os.path.basename(rules_path))[0]}_step{step_id}"

            ctl = clingo.Control(["0", "--opt-mode=optN", "--parallel-mode=36"])
            print(f"→ Traitement de {rules_name}")

            file_output = init_clingo(rules_path, rules_name)
            ctl.load(file_output)
            ctl.ground([("base", [])])
            ctl.solve(on_model=lambda model: on_model(model, rules_name))

        nb_total += nb_contigs
        step_id += 1

__main__()
