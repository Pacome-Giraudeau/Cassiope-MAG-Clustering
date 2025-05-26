

% Chaque contig est associé à un seul cluster
{ assigned(P, C) : cluster(C) } = 1 :- contig(P).


% Cohérence taxonomique : deux contigs d'un même cluster doivent partager un taxon au même rang
same_taxo(C1, C2, Rank) :- assigned(C1, Cluster), assigned(C2, Cluster), 
                            taxo(C1, Rank, T), taxo(C2, Rank, T), C1 != C2.

% Maximiser la cohérence taxonomique dans les clusters
#maximize { Rank @ 2 : same_taxo(C1, C2, Rank) }.

% Récupérer les rangs par cluster directement à partir des affectations
cluster_rank(C, R) :- assigned(P, C), rank(P, R).

% On identifie les min et max correctement avec des variables sûres
not_min_rank(C, R) :- cluster_rank(C, R), cluster_rank(C, R2), R2 < R.
min_rank(C, R) :- cluster_rank(C, R), not not_min_rank(C, R).

not_max_rank(C, R) :- cluster_rank(C, R), cluster_rank(C, R2), R2 > R.
max_rank(C, R) :- cluster_rank(C, R), not not_max_rank(C, R).

% Contrainte : tous les rangs entre min et max doivent exister
:- min_rank(C, Min), max_rank(C, Max), R = Min..Max, not cluster_rank(C, R).


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




#show assigned/2.