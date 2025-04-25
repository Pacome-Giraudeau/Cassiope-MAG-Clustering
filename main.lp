#script (python)

import clingo
import sys
import os
import numpy as np
from scipy.spatial.distance import euclidean

###### DESCRPITION DES CHEMINS VERS LES DIFF FICHIERS

folder = os.path.dirname(os.path.abspath(sys.argv[0]))
folder_data = "data"
folder_data_set = "first_set"

tokens_files = [
    os.path.join(folder, folder_data, folder_data_set, "scg.expected")
]

contigs_tokens_files = [
    os.path.join(folder, folder_data, folder_data_set, "contig_scg_details.tsv")
]

contigs_kmere_files = [
    os.path.join(folder, folder_data, folder_data_set, "kmer_contigs.clean_very_short.txt")
]

contigs = list()
tokens = dict()
contigs_dist = dict()
contigs_kmere = dict()
contigs_tokens = dict()


nb_clusters = 3

####### LECTURE DES DONN2ES

def read_tokens():
    for t_file in tokens_files:
        with open(t_file, "r", encoding="utf-8") as file:
            for line in file:
                try:
                    tokens[line.split()[0]].append(line.split()[1])
                except KeyError:
                    tokens[line.split()[0]] = [line.split()[1]]
                
def read_contig_scg():
    for scg_file in contigs_tokens_files:
        with open(scg_file, "r",  encoding="utf-8") as file:
            for line in file:
                contigs_tokens[line.split()[0]] = line.split()[1:12]
                print(line.split()[0], " : ", line.split()[1:12])

def read_contig_kmere():
    for kmere_file in contigs_kmere_files:
        with open(kmere_file, "r",  encoding="utf-8") as file:
            for line in file:
                if(line.split()[0] != "contig"):
                    contigs.append(line.split()[0])
                    contigs_kmere[line.split()[0]] = np.array(line.split()[1:])
                    contigs_kmere[line.split()[0]] = contigs_kmere[line.split()[0]].astype(int)

def calculate_distance_contigs():
    c = 0
    for c1 in contigs:
        c += 1
        print(c,"/",len(contigs))
        for c2 in contigs:
            contigs_dist[(c1, c2)] = euclidean(contigs_kmere[c1], contigs_kmere[c2])



######### CLINGO 

# Create a Clingo control object

def add_clingo_contigs(prg):
    for c1 in contigs:
        prg.add("base", [], f'contig("{c1}").')
        

def add_clingo_distance(prg):
    for (c1, c2) in contigs_dist.keys():
        if c1<c2:
            d = int(contigs_dist[(c1, c2)])
            prg.add("base", [], f'distance("{c1}","{c2}", {d}).')

def add_clingo_scg(prg):
    for t in tokens.items():
        prg.add("base", [], f'token("{t}").')
    for ct in contigs_tokens.keys():
        prg.add("base", [], f'contigs_Token("{ct}", {contigs_tokens[ct][0]}).')        

def add_clingo_clusters(prg):
    for cluster in range(nb_clusters):
        prg.add("base", [], f'cluster({cluster}).')     
        

def add_clingo_rules_pacome():
    return
            # distance_between_clusters(C1, C2, D) :-
            #     assigned(P1, C1), assigned(P2, C2),
            #     C1 != C2,
            #     distance(P1, P2, D).
            # min_distance(C1, C2, MinD) :-
            #     cluster(C1), cluster(C2), C1 < C2,
            #     MinD = #min { D : distance_between_clusters(C1, C2, D) }.
            
    
def init_clingo(prg):
    add_clingo_clusters(prg)
    add_clingo_contigs(prg)
    add_clingo_distance(prg)
    add_clingo_scg(prg)
    add_clingo_rules_pacome()


####### MAIN

def main(prg):
    
    read_contig_kmere()
    calculate_distance_contigs()
    
    # Ground and solve
    
    
    init_clingo(prg)
    prg.ground([("base", [])])
    prg.solve()
    
    return
    
#end.

{ assigned(P, C) : cluster(C) } = 1 :- contig(P).
total_dist(Dtot) :- 
    Dtot = #min { D : 
       assigned(P1, C), 
       assigned(P2, C), 
       distance(P1, P2, D), 
       cluster(C),
       P1 < P2}.
#minimize { D : 
   assigned(P1, C), 
   assigned(P2, C), 
   distance(P1, P2, D), 
   cluster(C),
   P1 < P2}.

% Cohérence taxonomique : deux contigs d'un même cluster doivent partager un taxon au même rang
% same_taxo(C1, C2, Rank) est vrai si C1 et C2 sont dans le même cluster et partagent un taxon au rang Rank
same_taxo(C1, C2, Rank) :- assigned(C1, Cluster), assigned(C2, Cluster), taxo(C1, Rank, T), taxo(C2, Rank, T), C1 != C2.

% Interdiction des clusters incohérents : si deux contigs ont des taxonomies différentes aux niveaux 1, 2 ou 3, ils ne peuvent pas être dans le même cluster
% a ajuster
:- assigned(C1, Cluster), assigned(C2, Cluster), taxo(C1, Rank, T1), taxo(C2, Rank, T2), C1 != C2, T1 != T2, Rank <= 1.

% Maximiser la cohérence taxonomique dans les clusters
% Maximiser la valeur de rank dans les clusters 
#maximize { Rank : same_taxo(C1, C2, Rank) }.


#show assigned/2.
#show contig/1.
#show taxo/3.
#show cluster/1.