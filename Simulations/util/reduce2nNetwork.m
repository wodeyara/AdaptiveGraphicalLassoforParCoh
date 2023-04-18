function newG1 = reduce2nNetwork(newG)
% given an augmented binarized network newG, this function reduces the
% network down to the normal cross-spectral density size, but still
% binarized. Note that if there are edges in the complementary matrix,
% this function will include them in the final binarized matrix. 
% ==============================
% INPUT
% newG: 2n x 2n matrix where n is number of nodes
% ==============================
% OUTPUT
% newG1: n x n matrix where n is number of nodes

newG1 = (newG(1:length(newG)/2 , 1:length(newG)/2)  + ...
         newG(1:length(newG)/2 , 1 + (length(newG)/2) : length(newG)) +...
            newG(1 + (length(newG)/2) : length(newG) , 1:length(newG)/2) + ... 
            newG(1 + (length(newG)/2) : length(newG) , 1 + (length(newG)/2) : length(newG))) > 0;