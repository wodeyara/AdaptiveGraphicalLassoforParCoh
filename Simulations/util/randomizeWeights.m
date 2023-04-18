function G = randomizeWeights(G1)
% this function takes any weighted graph G1 and shuffles the weights
% while keeping the same adjacency structure to generate the random weighted graph G 
% ===========================================
% INPUT
% G1: weighted graph 
% ===========================================
% OUTPUT
% G: randomized weighted graph 
%
% ani wodeyar

G1 = G1.*(~double(eye(length(G1)))); % getting rid of diagonal

if issymmetric(G1)
    locs = find(triu(abs(G1)>0,1));
    weights = G1(locs);
    r = randperm(length(locs));
    G = zeros(length(G1));
    G(locs) = weights(r);
    G = G+ G';
else
    locs = find(abs(G1)>0);
    weights = G1(locs);
    r = randperm(length(locs));
    G = zeros(length(G1));
    G(locs) = weights(r);
end