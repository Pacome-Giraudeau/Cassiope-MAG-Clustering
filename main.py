import clingo
import sys
import os
import numpy as np
from scipy.spatial.distance import euclidean

###### DESCRPITION DES CHEMINS VERS LES DIFF FICHIERS

folder = os.path.dirname(os.path.abspath(sys.argv[0]))
folder_data = "data"
folder_data_set = "first_set"
folder_data_set_test = "test_set"
folder_data_working = os.path.join(folder, folder_data, folder_data_set)

""" Description du chemin vers scg.expected qui contient 
la liste des gènes marqueurs bactériens attendus
dans un génome"""
tokens = dict() # Association gène marqueur, groupe taxonomique
tokens_files = [
    os.path.join(folder, folder_data, folder_data_set, "scg.expected")
]

""" Description du chemin vers mean_coverage_Q2Q3_contigs_2 qui contient 
les informations d'abondance
par échantillon de chacun des contigs """
contigs_coverage = dict()
contigs_coverage_files = os.path.join(folder, folder_data, folder_data_set, "mean_coverage_Q2Q3_contigs_3.txt")

""" Description du chemin vers contig_scg_details qui contient 
les gènes marqueurs identifiés pour
chaque contig et leurs affiliations taxonomiques """
contigs_tokens = dict() # Association contigs, gène marqeurs, groupe taxonomique détaillé
contigs_tokens_files = [
    os.path.join(folder, folder_data, folder_data_set, "contig_scg_details.tsv")
]

""" Description du chemin vers kmer_contigs.clean_very_short qui contient 
une version raccourcie de la table des fréquences en
tétranucléotides de chacun des contigs """
contigs_kmere = dict() # Décompte des bases tétranucléotide de chaque contig
contigs_kmere_files = [
    os.path.join(folder, folder_data, folder_data_set, "kmer_contigs.clean_short.txt")
]

contigs = list() # Initialisation de la liste des contigs


nb_clusters = 5

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
    n = len(contigs)
    
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

    ponderation = 1

    # Chargement des données
    contigs_coverage, liste_coverage = read_contig_coverage()
    contigs_kmere, liste_kmere = read_contig_kmere()

    # Identification des contigs en commun
    contigs_en_commun = list(set(liste_coverage) & set(liste_kmere))
    print("Liste des contigs en commun :", contigs_en_commun)
    for contig in contigs_en_commun :
        contigs.append(contig)
    n = len(contigs)
    
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

def get_distance_contigs(c1, c2):
    """
    c1 et c2 deux entiers
    Renvoie la distance entre c1 et c2
    """
    return int(contigs_dist[min(c1, c2), max(c1, c2)])
######### CLINGO 


def add_clingo_contigs():
    """Ajoute les contigs à Clingo.
    Input : Liste des contigs
    Output : Règles Clingo contig("contig1")."""

    print("Vérification de la liste des contigs :", contigs)

    for c1 in contigs:
        print(f"Ajout du contig {c1} dans Clingo")
        ctl.add("base", [], f'contig("{c1}").')

        

def add_clingo_distance():
    """Ajoute les distances entre les contigs à Clingo.
    Input : Matrice des distances, Dictionnaire {contig : index}
    Output : Règles Clingo distance("contig1", "contig2", 5)."""
    print("Ajout des distances entre contigs...")
    for c1 in range(len(contigs)):
        
        print("[Ajout des distances entre contigs [", c1, "/", len(contigs)-1,"]")
        for c2 in range(c1, len(contigs)):
            ctl.add("base", [], f'distance("{contigs[c1]}","{contigs[c2]}", {get_distance_contigs(c1, c2)}).')
# Sortie 
# distance("contig1", "contig2", 5).

def add_clingo_scg():
    """Ajoute les gènes marqueurs et la taxonomie des contigs à Clingo.
    Input : Dictionnaire {contig : ([gènes marqueurs], [taxonomie])}
    Output : Règles Clingo scg("contig1", "gene1")."""

    '''for t in tokens.items():
        ctl.add("base", [], f'token("{t}").')
    #for ct in contigs_tokens.keys():
        #ctl.add("base", [], f'contigs_Token("{ct}", {contigs_tokens[ct][0]}).') '''
        
    print("Ajout des contigs marqueurs et de leurs taxonomie...")
    for t in tokens.items():
        ctl.add("base", [], f'token("{t}").')
    
    for entry in contigs_tokens:
        contig = entry[0]
        genes = entry[1]
        taxonomie = entry[2]

    for gene in genes:
        ctl.add("base", [], f'scg("{contig}", "{gene}").')
        print("Ajout SCG Clingo :", f'scg("{contig}", "{gene}").')

    if taxonomie:
        taxo_string = ", ".join(f'"{level}"' for level in taxonomie)
        ctl.add("base", [], f'taxo("{contig}", {taxo_string}).')
        print("Ajout Taxonomie Clingo :", f'taxo("{contig}", {taxo_string}).') 

def add_clingo_scg_and_taxonomy():
    """Ajoute les gènes marqueurs et la taxonomie des contigs à Clingo.
    Input : Dictionnaire {contig : ([gènes marqueurs], [taxonomie])}
    Output : Règles Clingo scg("contig1", "gene1")."""
    

    print("Ajout des contigs marqueurs et de leurs taxonomie...")
    for contig, (genes, taxo) in contigs_tokens.items():
        for gene in genes:
            #print(f'Ajout Clingo : scg("{contig}", "{gene}")')
            ctl.add("base", [], f'scg("{contig}", "{gene}").')

        for rank, taxon in enumerate(taxo, start=1):
            #print(f'Ajout Clingo : taxo("{contig}", {rank}, "{taxon}")')
            ctl.add("base", [], f'taxo("{contig}", {rank}, "{taxon}").')

def add_clingo_clusters():
    """Ajoute les clusters à Clingo.
    Input : Nombre de clusters  
    Output : Règles Clingo cluster(1)."""

    for cluster in range(nb_clusters):
        ctl.add("base", [], f'cluster({cluster}).') 

        
def add_clingo_rules():
    
    print("Ajout des règles  basiques...")
    ctl.add("base", [],"""
        % Chaque contig est associé à un seul cluster
        { assigned(P, C) : cluster(C) } = 1 :- contig(P).
        
        """)
    
def add_clingo_rules_pacome():
    print("Ajout des règles de Pacôme...")
    ctl.add("base", [],"""
            distance_clusters(Cluster1, Cluster2, Dc) :- 
                assigned(C1, Cluster1), assigned(C2, Cluster2), 
                cluster(Cluster1), 
                cluster(Cluster2), 
                Cluster1 < Cluster2, 
                distance(C1, C2, Dc).



            min_distance(Cluster1, Cluster2, Dc1) :- 
                distance_clusters(Cluster1, Cluster2, Dc),
                distance_clusters(Cluster1, Cluster2, Dc1), 
                cluster(Cluster1), 
                cluster(Cluster2), 
                Cluster1 < Cluster2,
                Dc1 <= Dc.
                
            #maximize { Dc@2 : min_distance(Cluster1, Cluster2, Dc)}.

            #minimize { D@2 : 
                assigned(C1, Cluster), 
                assigned(C2, Cluster), 
                distance(C1, C2, D), 
                cluster(Cluster),
                C1 < C2}.
            """)
# Règles :  chaque contig est associé à seulement un cluster
#           la somme totale des distances associées aux contigs doit être minimisée
#           comptage de la complétion

def add_clingo_rules_sarra():
    print("Ajout des règles de Sarra...")
    '''ctl.add("base", [], """
        % Cohérence taxonomique : deux contigs d'un même cluster doivent partager un taxon au même rang
        same_taxo(C1, C2, Rank) :- assigned(C1, Cluster), assigned(C2, Cluster), taxo(C1, Rank, T), taxo(C2, Rank, T), C1 != C2.

        % Interdiction des clusters incohérents : si deux contigs ont des taxonomies différentes aux niveaux 1, 2 ou 3, ils ne peuvent pas être dans le même cluster
        :- assigned(C1, Cluster), assigned(C2, Cluster), taxo(C1, Rank, T1), taxo(C2, Rank, T2), C1 != C2, T1 != T2, Rank <= 3.

        % Maximiser la cohérence taxonomique dans les clusters
        #maximize { Rank : same_taxo(C1, C2, Rank) }.

        % Afficher uniquement les clusters et leurs contenus
        #show assigned/2.
    """)'''
    ctl.add("base", [], """
        % Maximiser la complétion 
        gene_in_cluster(G, B) :- assigned(C, B), scg(C, G).

        % Calcul du nombre de gènes distincts dans chaque cluster
        count_genes(B, Count) :- cluster(B), Count = #count { G : gene_in_cluster(G, B) }.

        % On maximise le total des gènes distincts sur chaque cluster
        #maximize { Count, B : count_genes(B, Count) }.

        % Minimiser la redondance 
        gene_frequency(G, B, Freq) :- assigned(_, B), scg(_, G), cluster(B),
                                       Freq = #count { C : assigned(C, B), scg(C, G) }.

        #minimize { Freq - 1, G, B : gene_frequency(G, B, Freq), Freq > 1 }.

    """)


def add_clingo_rules_marion():
    print("Ajout des règles de Marion...")
    ctl.add("base", [], """
               
        % --- Ajout du calcul du diamètre maximum d'un cluster ---
        diam(B, D) :- assigned(C1, B), assigned(C2, B), C1 < C2, cluster(B), distance(C1, C2, D).
        diam(B, D-1) :- diam(B, D), D > 0.
        maxDiameter(B, D) :- diam(B, D), not diam(B, D+1).


        #minimize { D@2 : maxDiameter(B, D) }.

        % Cohérence taxonomique : deux contigs d'un même cluster doivent partager un taxon au même rang
        same_taxo(C1, C2, Rank) :- assigned(C1, Cluster), assigned(C2, Cluster), 
                                   taxo(C1, Rank, T), taxo(C2, Rank, T), C1 != C2.

        % Interdiction des clusters incohérents : si deux contigs ont des taxonomies différentes aux niveaux 1, 2 ou 3,
        % ne peuvent pas être dans le même cluster
        :- assigned(C1, Cluster), assigned(C2, Cluster), 
           taxo(C1, Rank, T1), taxo(C2, Rank, T2), 
           C1 != C2, T1 != T2, Rank <= 3.

        % Maximiser la cohérence taxonomique dans les clusters
        #maximize { Rank@1 : same_taxo(C1, C2, Rank) }.

    """)


'''def on_model(model):
    
    if(model.optimality_proven):
        """Callback function that is called for each solution."""
    atoms = model.symbols(shown=True)
    filtered_atoms = [str(atom) for atom in atoms if atom.name == "assigned" or atom.name == "total_dist" ] # atom.name == "assigned" or 
    
    print("model : ", model.number, " ( opti : ", model.optimality_proven , " )  Solution optimale:")
    print(" \n".join(filtered_atoms))'''


def on_model(model):
    global model_count
    model_count += 1

    atoms = model.symbols(shown=True)
    count_atoms = 0
    for atom in atoms:
        count_atoms += 1

    # Récupérer uniquement les affectations aux clusters et le nombre de contigs par cluster
    assigned_atoms = [str(atom) for atom in atoms if atom.name == "assigned"]
    count_assigned_atoms = len(assigned_atoms)
    cluster_count_atoms = []
    clusters = [int(atom.split(',')[-1][:-1]) for atom in assigned_atoms]
    # print(clusters)
    for cluster in clusters:
        if cluster not in cluster_count_atoms:
            cluster_count_atoms.append(cluster)



    # Compter le nombre de clusters
    cluster_count = len(cluster_count_atoms)

    print(f"Nombre de clusters : {cluster_count}")
    print(f"Nombre de contigs total : {count_atoms}")
    print(f"Nombre de contigs assignés : {count_assigned_atoms}")
    print("\nRésultat du clustering :")
    print("\n".join(assigned_atoms) if assigned_atoms else "Aucun contig assigné")


    print("\nFin de l'affichage des clusters\n")
    print("\n optimal ? ", model.optimality_proven)

    with open(os.path.join(folder_data_working, "clusters.txt"), "w") as f:
        for atom in atoms:
            if atom.name == "assigned" and len(atom.arguments) == 2:
                contig, cluster = atom.arguments
                contig = str(contig).strip('"')
                f.write(f"{contig}\tbin_{cluster}\n")

    if model_count >= 5:
        return False


def init_clingo():
    add_clingo_clusters()
    add_clingo_contigs()
    add_clingo_distance()
    add_clingo_scg_and_taxonomy()
    add_clingo_rules()
    add_clingo_rules_marion()
    add_clingo_rules_pacome()
    add_clingo_rules_sarra()


####### MAIN

contigs_dist, contigs = calculate_distance_contigs()
    
# Create a Clingo control object
ctl = clingo.Control(["0", "--opt-mode=opt"]) # "--parallel-mode=12"
def __main__():

    global model_count
    model_count = 0
    
    read_contig_kmere()
    read_contig_scg()
    
    # Ground and solve
    print("\nFin de la phase d'initialisation\n")
    
    
    ctl.configuration.verbose = 2
    ctl.configuration.opt_mode = "optN"
    
    init_clingo()

    print("Initialisation de clingo terminée !")
    ctl.ground([("base", [])]) 
    print("Grounding de clingo terminée !")
    ctl.solve(on_model=on_model)
    print("Résolution de clingo terminée !")
    
    return
    
    
__main__()



