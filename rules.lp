

% Chaque contig est associé à un seul cluster
{ assigned(P, C) : cluster(C) } = 1 :- contig(P).


% Cohérence taxonomique : deux contigs d'un même cluster doivent partager un taxon au même rang
same_taxo(C1, C2, Rank) :- assigned(C1, Cluster), assigned(C2, Cluster), 
                            taxo(C1, Rank, T), taxo(C2, Rank, T), C1 != C2.

% Interdiction des clusters incohérents : si deux contigs ont des taxonomies différentes aux niveaux 1, 2 ou 3,
% ne peuvent pas être dans le même cluster
:- assigned(C1, Cluster), assigned(C2, Cluster), 
    taxo(C1, Rank, T1), taxo(C2, Rank, T2), 
    C1 != C2, T1 != T2, Rank <= 3.

% Maximiser la cohérence taxonomique dans les clusters
#maximize { Rank : same_taxo(C1, C2, Rank) }.


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
    

% Maximisation de la somme des distances entre BINS
#maximize { Dc : min_distance(Cluster1, Cluster2, Dc)}.

% Minimisation somme des distances intra
#minimize { D : 
    assigned(C1, Cluster), 
    assigned(C2, Cluster), 
    distance(C1, C2, D), 
    cluster(Cluster),
    C1 < C2}.




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



% --- Ajout du calcul du diamètre maximum d'un cluster ---
diam(B, D) :- assigned(C1, B), assigned(C2, B), C1 < C2, cluster(B), distance(C1, C2, D).
diam(B, D-1) :- diam(B, D), D > 0.
maxDiameter(B, D) :- diam(B, D), not diam(B, D+1).


#minimize { D : maxDiameter(B, D) }.


#show assigned/2.