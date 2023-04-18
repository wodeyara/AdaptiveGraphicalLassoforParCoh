function [allS, weightingMatrix]  = getDipOrientations(L, data, regParam, oldWeightingMat)
% this function uses the data to define the dipole orientations given that
% the lead field contains the lead field for the x/y/z directiosn at a
% number of ROIs. 
% L = channels x (sources*3)
% data = channels x time
% regParam is the percentile used for inversemodel (line 13)
% Ani Wodeyar

    if nargin<4
        oldWeightingMatFlag =0;
        oldWeightingMat =[];
    else
        oldWeightingMatFlag =1;
    end
    W = diag(1./sum((L.^2))); %||l_i||^2 summing over channels
    W = (W); % based on Hamalainens work from 2006 Neuroimage

    opInv = inversemodel(L*W, 'prctile', regParam);

    dataRecon = opInv * data;
    weightingMatrix = zeros(size(L,2)/3, size(L,2));
    allS = zeros(size(L,2)/3,1);
    dataRecon = dataRecon';
    for i = 3:3:size(L,2)
        tmpDipole = dataRecon(:,i-2:i);
        [~,S,V] = svd(tmpDipole, 'econ');
        S = diag(S.^2)/sum(diag(S.^2));
        if oldWeightingMatFlag
            oldV = oldWeightingMat(i/3, i-2:i);
            if round(dot(oldV, V(:,1))) == -1 % want to keep sign consistent across multiple weighting matrices
                V(:,1) = -V(:,1);
            end            
        end
        weightingMatrix(i/3, i-2:i) = V(:,1);
        allS(i/3) = S(1);
    end