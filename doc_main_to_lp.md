# Documentation de main_to_lp

## Descriion générale

main_to_lp permet de clusteriser des contigs via le langage de programmation ASP et son API python

## Requis
### Imports
Clingo

'''Inserer tuto'''

### Bibilothèques pythons
clingo
sys
os
numpy
scipy
subprocess
'''Inserer tuto'''



## Format des données

- scg.expected  

Description du chemin vers scg.expected qui contient la liste des gènes marqueurs bactériens attendus dans un génome.

- mean_coverage_Q2Q3_contigs.txt

Description du chemin vers mean_coverage_Q2Q3_contigs qui contient les informations d'abondance par échantillon de chacun des contigs.

- contig_scg_details.tsv

Description du chemin vers contig_scg_details qui contient les gènes marqueurs identifiés pour chaque contig et leurs affiliations taxonomiques.

- kmer_contigs.clean_short

Description du chemin vers kmer_contigs.clean_very_short qui contient une version raccourcie de la table des fréquences entétranucléotides de chacun des contigs.

## Format des données de sortie
XX est le nombre de contigs, c'est là pour différencier les tests.

#### XX_clusters.lp
Ensembles des faits liés aux clusters. 
Format : 
    
    cluster(C). #où C est un entier, ie le numéro de cluster

#### XX_distances.lp
Ensembles des faits liés aux distances entre contigs. 
Format : 
    
    ditance("contig1", "contig1", D). #où D est un entier

#### XX_scg_&_taxonomy.lp
Ensembles des faits liés aux taxonomies et aux gènes marqueurs. 
Format : 
    
    cscg("contig1", "gene1").
    taxo("{contig}", {rank}, "{taxon}").

#### XX_contigs.lp
Ensembles des faits liés aux contigs. 
Format : 
    
    contig("contig"). 

#### rules.lp
Ensembles des règles d'optimisation. 

## Paramètres

nb_clusters : le nombre de clusters*

config_parallele = "--parallel-mode=" + str(10)
config_opti = "--opt-mode=opt"
config_verbose = "--verbose=0"

## Fonctions 

##### LECTURE DES DONN2ES

#### def read_tokens():
    
    """ 
    Lecture du fichier scg.expected et extraction des gènes marqueurs bactériens.
    Input : Fichier scg.expected
    Output : Dictionnaire {gène marqueur : groupe taxonomique}
    """
                

#### def read_contig_scg():
    
    """
    Lecture du fichier contig_scg_details.tsv et extraction des données taxonomiques + gènes marqueurs.
    Input : Fichier contig_scg_details.tsv
    Output : Dictionnaire {contig : ([gènes marqueurs], [taxonomie])}
    """

#### def read_contig_kmere():
    
    """
    Lecture du fichier kmer_contigs.clean.txt et extraction des données tétranucléotides.
    Input : Fichier kmer_contigs.clean
    Output : Dictionnaire {contig : [fréquences tétranucléotides]}, Liste [contigs]"""


#### def read_contig_coverage():
    
    """
    Lecture du fichier mean_coverage_Q2Q3_contigs_2.txt et extraction des données d'abondance.
    Input : Fichier mean_coverage_Q2Q3_contigs_2.txt
    Output : Matrice (contigs x stations), Liste [contigs]
    """


#### def calculate_distance_contigs_naif():
    """
    Calcul de la distance euclidienne entre les contigs à partir des données d'abondance et de tétranucléotides.
    Création d'un vecteur combiné pour chaque contig et calcul de la distance euclidienne entre chaque paire de contigs.
    Input : contigs_coverage, contigs_kmere
    Output : Matrice des distances, Dictionnaire {contig : index}
    """


#### def calculate_distance_contigs():

    """
    Calcul de la distance euclidienne entre les contigs à partir des données d'abondance et de tétranucléotides.
    Calcul de la distance euclidienne entre les vecteurs d'abondance puis de tétranucléotides de chaque paire de contigs
    puis moyenne avec pondération.
    Input : contigs_coverage, contigs_kmere
    Output : Matrice des distances, Dictionnaire {contig : index}
    """


#### def get_distance_contigs(c1, c2, contigs_dist):

    """
    c1 et c2 deux entiers
    Renvoie la distance entre c1 et c2
    """


#### def get_nb_contigs():
    
    """
    renvoie le nombre de contigs considérés
    """


#### def get_contigs():
    
    """
    renvoie la liste des contigs
    """
    

##
##### CLINGO 
##


#### def add_clingo_contigs():
    
    """
    Ajoute les contigs à Clingo.
    Input : Liste des contigs
    Effet de bord : Règles Clingo contig("contig1").
    - Si n'existe pas création du fichier contigs.lp avec les règles
    Output : chemin vers le fichier mentionné ci dessus
    """

    
#### def add_clingo_distance():
    
    """
    Ajoute les distances entre les contigs à Clingo.
    Input : Matrice des distances, Dictionnaire {contig : index}
    Effet de bord : Règles Clingo distance("contig1", "contig2", 5).
    - Si n'existe pas création du fichier distances.lp avec les distances
    Output : chemin vers le fichier mentionné ci dessus
    """
    

#### def add_clingo_scg_and_taxonomy():
    
    """
    Ajoute les gènes marqueurs et la taxonomie des contigs à Clingo.
    Input : Dictionnaire {contig : ([gènes marqueurs], [taxonomie])}
    Effet de bord : Règles Clingo 
              - scg("contig1", "gene1").
              - taxo("{contig}", {rank}, "{taxon}").
    
    - Si n'existe pas création du fichier token_&_taxo.lp avec les tokens et les taxonomoies
    Output : chemin vers le fichier mentionné ci dessus
    """

#### def add_clingo_clusters():
    
    """
    Ajoute les clusters à Clingo.
    Input : Nombre de clusters  
    Output : Règles Clingo cluster(1).
    - Si n'existe pas création du fichier clusters.lp avec les clusters
    - load des clusters dans le controleur clingo
    """


#### def add_clingo_rules():
    
    """
    Ajoute les clusters à Clingo.
    Input : Nombre de clusters  
    Output : Règles Clingo cluster(1).
    - Dans tous les cas création du fichier rules.lp avec les règles
    - load des règles dans le controleur clingo
    """



#### def on_model(model):
    
    """Fonction d'affichage des models"""


#### def init_clingo():

    """
    Cration des 5 fichiers 
    - clusters.lp
    - contigs.lp
    - distances.lp
    - scg_&_taxo.lp
    - rules.lp
   
    Merging de ces 5 fichiers programme en 1 fichier programme unique
    - output.lp

    Output : chemin vers le fichier programme output.lp
    
    """




