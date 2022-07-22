%% Getting the connectivity and covariance: 
% here we get first the randomized weights precision structure 
% then we generate the covariance matrix assuming a multivariate normal process (here
% a complex valued multivariate normal process) and we get the covariance matrix on

useAreas = setdiff(unique(labels), subcortROIs);

phaseFromLength = (pi/2+randn(length(useAreas))*.25) * double(~eye(length(useAreas))) .*(SC(useAreas, useAreas)>0) ;
% first we will normalize each row sum to 1.
cnt  = 1;
err = 1;
while err ~=0 && cnt < 1e3
    r = abs(100+randn(720,1)*30);
    G = double(triu(SC(useAreas, useAreas)>0,1));
    G(G>0) = r;
    G = G + G';
    orig_G = G;
    origNetwork = abs(orig_G)>0;
    G = triu(G .*exp(1i * phaseFromLength),1); % phase from length comes from the fiberTrackLengths.m script currently in lesion_FC_SC.
    G = G+G';
    Q = real2Complex([G,zeros(length(G))*sqrt(-1); zeros(length(G))*sqrt(-1), conj(G)], 1); % getting real valued G
    
    tmpDegs = sum(abs(Q.*double(~eye(228))),2);
    [Q] = makePosDef(diag(tmpDegs) + Q); % this step is akin to solving a linear model
    % now get the covariance structure for this process.
    Q_compVal = real2Complex(Q,0);
    Q_compVal = Q_compVal(1:length(useAreas), 1:length(useAreas));

    [~,err] = chol(Q);
    cnt  = cnt +1;
end
covMat = Q\eye(size(Q));% helps with numerical error
data = mvnrnd(zeros(length(covMat),1), covMat,samps)';
compData = data(1:length(useAreas),:) +1j*data(length(useAreas)+1:2*length(useAreas),:);

covMatCSD = real2Complex(covMat, 0);
covMatCSD = covMatCSD(1:length(useAreas), 1:length(useAreas));

orig_ParCohMat = abs(Q_compVal ./sqrt(diag(Q_compVal) * diag(Q_compVal)'));

clear data compData orig_CohMat orig_CovMat r frequency axonalSpeed phaseFromLength cnt err covMatCSD
