% Chaque contig est associé à un seul cluster
{ assigned(P, C) : cluster(C) } = 1 :- contig(P).


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


#show assigned/2.
