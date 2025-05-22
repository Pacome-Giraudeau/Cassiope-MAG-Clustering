import clingo
import sys
import os
import numpy as np
from scipy.spatial.distance import euclidean
import subprocess
import random
import glob
from Bio import Phylo
from io import StringIO

###### DESCRPITION DES CHEMINS VERS LES DIFF FICHIERS

nb_clusters =4
nb_contigs_to_classify=200

folder = os.path.dirname(os.path.abspath(sys.argv[0]))
folder_data = "data"
folder_data_set = "first_set"
folder_data_set_test = "test_set"
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

copy_random_lines_with_scg(folder_data_working=folder_data_working, nb_contigs=nb_contigs_to_classify)

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
contigs_coverage_files = os.path.join(folder, folder_data, folder_data_set, f"{nb_contigs_to_classify}_mean_coverage_Q2Q3_subset.txt")

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
    os.path.join(folder, folder_data, folder_data_set, f"{nb_contigs_to_classify}_kmer_contigs_subset.txt")
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

    ponderation = 100
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


    # Calcul de la matrice des distances
    for i in range(n):
        c1 = contigs[i]
        print(i + 1, "/", n)

        # Construction du vecteur combiné pour c1
        coverage_c1 = contigs_coverage[liste_coverage.index(c1)]
        kmere_c1 = contigs_kmere[c1]
        for j in range(i, n):  # Évite les doublons (distance symétrique)
            c2 = contigs[j]
            coverage_c2 = contigs_coverage[liste_coverage.index(c2)]
            kmere_c2 = contigs_kmere[c2]

            # Calcul de la distance euclidienne
            distance_coverage = euclidean(coverage_c1, coverage_c2)
            distance_kmere = euclidean(kmere_c1, kmere_c2)
            distance= (distance_coverage*ponderation + distance_kmere)/(ponderation+1)
            contigs_dist[i, j] = distance
            
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
    return nb - 1 # -1 car la première ligne, c'est "contig AAAA AACT ..."

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


# Parcours post-ordre de droite à gauche
def right_to_left_postorder(clade, ordered_contigs):
    if clade.is_terminal():
        ordered_contigs.append(clade.name)
    else:
        children = clade.clades
        if len(children) == 2:
            right_to_left_postorder(children[1], ordered_contigs)
            right_to_left_postorder(children[0], ordered_contigs)
        elif len(children) == 1:
            right_to_left_postorder(children[0], ordered_contigs)


def create_contigs_ranks():
    """Crée un fichier contig_ranks.tsv contenant les rangs des contigs à partir d'un arbre phylogénétique au format Newick.
    Input : Fichier ordre_item.txt (arbre phylogénétique)
    Output : Fichier contig_ranks.tsv (contig, rang)"""
    
    # Lecture du fichier ordre_item.txt
    with open(os.path.join(folder, "ordre_items.txt"), "r") as file:

        newick_str = file.read().strip()  # on retire les éventuels sauts de ligne inutiles

    # Parser l’arbre à partir de la chaîne Newick
    handle = StringIO(newick_str)
    tree = Phylo.read(handle, "newick")

    # Liste ordonnée des feuilles (contigs)
    ordered_contigs = []

    # Lancer le parcours
    right_to_left_postorder(tree.root, ordered_contigs)

    # Attribuer un rang (droite = rang 0)
    ranks = {contig: rank for rank, contig in enumerate(reversed(ordered_contigs))}

    # Sauvegarde dans un fichier tsv
    with open("ordre_contig.txt", "w") as f:
        for contig, rank in ranks.items():
            f.write(f"{contig}\t{rank}\n")

    print("Rangs des contigs enregistrés dans ordre_contig.txt.")

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



        

def add_clingo_distance():
    """Ajoute les distances entre les contigs à Clingo.
    Input : Matrice des distances, Dictionnaire {contig : index}
    Output : Règles Clingo distance("contig1", "contig2", 5)."""
    
    print("Création des distances...")
    
    file_path = os.path.join(folder_data_working, "programme_lp", str(nb_contigs) + "_distances.lp")
    if not os.path.exists(file_path):
        contigs_dist, contigs = calculate_distance_contigs()
        
        with open(file_path, "w") as f:
            for c1 in range(nb_contigs):
                
                print("[Ajout des distances entre contigs [", c1, "/", nb_contigs-1,"]")
                for c2 in range(nb_contigs):
                    f.write(f'distance("{contigs[c1]}","{contigs[c2]}", {get_distance_contigs(c1, c2, contigs_dist)}).\n')
    
    print("Création des distances effectué !\n")
    
    # ctl.load(file_path)
    return file_path


def add_clingo_ranks():
    """
    Ajoute les rangs des contigs à Clingo.
    Input : Fichier ordre_items.txt (tabulé : contig<TAB>rank)
    Output : Fichier LP contenant des règles Clingo de type rank("contig1", 5).
    """
    print("Création des rangs des contigs pour Clingo...")

    input_path = os.path.join(folder, "ordre_contig.txt")
    output_path = os.path.join(folder_data_working, "programme_lp", f"{nb_contigs}_ranks.lp")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    if not os.path.exists(output_path):
        with open(input_path, "r") as f_in, open(output_path, "w") as f_out:
            for line in f_in:
                line = line.strip()
                if not line or "\t" not in line:
                    continue  # ignore lignes vides ou malformées
                try:
                    name, rank = line.split("\t")
                    name = name.strip()
                    rank = int(rank.strip())
                    f_out.write(f'rank("{name}", {rank}).\n')
                except ValueError:
                    print(f"Ligne ignorée (mal formée) : {line}")

        print("Fichier de rangs Clingo créé :", output_path)
    else:
        print("Fichier de rangs Clingo déjà existant :", output_path)

    return output_path



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



def on_model(model, rules_name):
    global model_count
    model_count += 1

    atoms = model.symbols(shown=True)
    count_atoms = 0
    for atom in atoms:
        count_atoms += 1
    assigned_atoms = [str(atom) for atom in atoms if atom.name == "assigned"]
    count_assigned_atoms = len(assigned_atoms)
    cluster_count_atoms = []
    clusters = [int(atom.split(',')[-1][:-1]) for atom in assigned_atoms]
    for cluster in clusters:
        if cluster not in cluster_count_atoms:
            cluster_count_atoms.append(cluster)
    cluster_count = len(cluster_count_atoms)

    print(f"\nRésultat pour {rules_name} :")
    print(f"Nombre de clusters : {cluster_count}")
    print(f"Nombre de contigs total : {count_atoms}")
    print(f"Nombre de contigs assignés : {count_assigned_atoms}")
    print("\nRésultat du clustering :")
    print("\n".join(assigned_atoms) if assigned_atoms else "Aucun contig assigné")


    print("\nFin de l'affichage des clusters\n")
    print("\n optimal ? ", model.optimality_proven)

    # Sauvegarde
    clusters_output_file = os.path.join(folder_data_working, f"{nb_contigs_to_classify}_clusters_{rules_name}.txt")
    with open(clusters_output_file, "w") as f:
        for atom in atoms:
            if atom.name == "assigned" and len(atom.arguments) == 2:
                contig, cluster = atom.arguments
                contig = str(contig).strip('"')
                f.write(f"{contig}\tbin_{cluster}\n")

    if model_count >= 5:
        return False



def init_clingo(rules_file_path, rules_name):

    file_clusters = add_clingo_clusters()
    file_contigs = add_clingo_contigs()
    file_distances = add_clingo_distance()
    file_scg_taxo = add_clingo_scg_and_taxonomy()
    file_ranks = add_clingo_ranks()  # retourne maintenant un fichier .lp

    file_output = os.path.join(folder, f"output_{rules_name}.lp")
    print("Merging : ")
    print(file_clusters)
    print(file_contigs)
    print(file_distances)
    print(file_scg_taxo)
    print(rules_file_path)
    
    f1 = open(file_clusters, "r")
    f2 = open(file_contigs, "r")
    f3 = open(file_distances, "r")
    f4 = open(file_scg_taxo, "r")
    f5 = open(rules_file_path, "r")  
    f6 = open(file_ranks, "r")

    c1 = f1.read()
    c2 = f2.read()
    c3 = f3.read()
    c4 = f4.read()
    c5 = f5.read()
    c6 = f6.read()

    f7 = open(file_output, "w")

    f7.write(c1 + c2 + c3 + c4 + c5 + c6)

    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    f6.close()
    f7.close()


    print("Création du fichier progromme OK pour :", rules_name)

    return file_output



#LISTE DES FICHIERS RULES
rules_files = glob.glob(os.path.join(folder, "rules", "**", "rules_*.lp"), recursive=True)

print(f"Fichiers rules détectés :")
for rule in rules_files:
    print(" -", rule)

####### MAIN


nb_contigs = get_nb_contigs()
    
def __main__():

    

    global model_count
    model_count = 0

    
    read_contig_kmere()
    read_contig_scg()
    
    # Ground and solve
    print("\nFin de la phase d'initialisation\n")
    for rules_path in rules_files:
        model_count = 0
        
        rules_name = os.path.splitext(os.path.basename(rules_path))[0] 
        # Create a Clingo control object
        ctl = clingo.Control(["0", "--opt-mode=optN", "--parallel-mode=8"])

        print(f"\n==========\nTraitement pour {rules_name}\n==========\n")

        file_output = init_clingo(rules_path, rules_name)  
        
        print(f"Loading {file_output}...")
        ctl.load(file_output)
        print("Loading terminé !")
        print("Grounding...")
        ctl.ground([("base", [])]) 
        print("Grounding de clingo terminée !")
        print("Résolution...")
        result=ctl.solve(on_model=lambda model: on_model(model, rules_name))
        print("Resultat", result) 
        print("Résolution de clingo terminée !")
    
    return
    
    
__main__()



