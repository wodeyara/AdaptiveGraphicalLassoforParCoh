
% Currently only the weighted L2:
goodchans = 1:size(L,1);
W = 1./diag(sum((L.^2))); %||l_i||^2 summing over channels
W(isinf(W)) = 0;
W = (W.^powerForWeightedL2);

cnt =1;
allSourcePens = [1:91];
allOpInv = zeros(length(allSourcePens), size(L,1), size(L,2));
for penalization=   allSourcePens  
    [opInv] = inversemodel(L*W,'prctile',penalization );
    opInv = opInv';
    allOpInv(cnt , : , :) = opInv;
    cnt = cnt +1;
end

clear W opInv cnt allSourcePens