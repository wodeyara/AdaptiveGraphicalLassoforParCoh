clear all

% build a positive definite complex-valued precision matrix: 
S = [6, -.9*exp(1i*pi/4), -.5*exp(1i*pi/8),.55*exp(1i*pi/6), 0;
    -.9*exp(-1i*pi/8), 1.01, 0, 0, -.64*exp(1i*pi/8) ;
    -.5*exp(-1i*pi/8), 0, .6, 0, 0;
    .55*exp(-1i*pi/8), 0, 0, 1, 0;
    0,       -.64*exp(1i*pi/8), 0,0,1.5];
S = .5*(S+S');

% generate the augmented complex valued precision
S_aug = [S,zeros(5);zeros(5),conj(S)];
% convert to the real-valued complex precision
realPrec = real2Complex(S_aug,1); % flag 1 means the augmented matrix is complex-valued

% write out binarized prior matrix
origNetwork = [1,1,1,1,0;
               1,1,0,0,1;
               1,0,1,0,0;
               1,0,0,1,0;
               0,1,0,0,1];
origNetwork = origNetwork - eye(5);

% build a false network while ensuring that number of edges is the same
[fakeNetwork] = randomizer_bin_und(origNetwork,1);
numEdgesToAdd = sum(sum(triu(fakeNetwork.*origNetwork,1))); 
fakeNetwork = (fakeNetwork - fakeNetwork.*origNetwork); % remove SC edges if any
tmp1 = double(~fakeNetwork) - eye(length(fakeNetwork)) - origNetwork;
tmp1 = triu(tmp1,1);
tmp1Locs = find(tmp1);
tmp1LocsToUse = randperm(length(tmp1Locs));
tmp1LocsToUse = tmp1Locs(tmp1LocsToUse(1:numEdgesToAdd));
fakeNetwork(tmp1LocsToUse) = 1;
fakeNetwork = double((triu(fakeNetwork,1) + triu(fakeNetwork,1)') > 0);
fakeNetwork = fakeNetwork + eye(length(fakeNetwork));

% network for the graphical lasso:
GforFit =[(eye(length(origNetwork)) +double(origNetwork)), ... 
    double(origNetwork) ; double(origNetwork), (eye(length(origNetwork)) +double(origNetwork))];
     
% penalization parameters
lambdasIn = fliplr([.9,.8,.7,.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001]);           
allLambdasOut = fliplr([.9,.8,.7,.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001]); %           
samps = 240;
% SNR = 10;
%%
parfor i = 1:200
    tic
    %simulate data from precision:
    data = mvnrnd(zeros(length(realPrec),1),(realPrec)\eye(size(realPrec)),samps)';
    
    % separate data into ensembles
    reconDataReal = permute(reshape(data', 4,samps/4, size(data,1)),[1,3,2]);

    % estimate precision using adaptive graphical lasso with correct network
    [networkPrecComp, penInComp(i), penOutComp(i)] = estBestPenalizationQUIC(... 
        reconDataReal, origNetwork,lambdasIn,allLambdasOut, 0);
    
    % estimate precision using adaptive graphical lasso with fake network
    [networkPrecCompFake, penInCompFake(i), penOutCompFake(i)] = estBestPenalizationQUIC(... 
        reconDataReal, fakeNetwork,lambdasIn,allLambdasOut, 0);    
    
    % check recovery under each approach above:
    corrsNet(i) =(corr((realPrec(triu(GforFit>0,1)|triu(networkPrecComp>0,1)) + 1e-4) ...
                     ,(networkPrecComp(triu(GforFit>0,1)|triu(networkPrecComp>0,1))+ 1e-4)));
    newG1 = reduce2nNetwork(abs(networkPrecComp)>0);
                 
    edgesInNetwork(i) = sum(sum(newG1.*triu(origNetwork,1)));
    edgesNotInNetwork(i) =  sum(sum(newG1.*triu(~origNetwork,1)));
 
    corrsNetFake(i) =(corr((realPrec(triu(GforFit>0,1)|triu(networkPrecCompFake>0,1)) + 1e-4) ...
                     ,(networkPrecCompFake(triu(GforFit>0,1)|triu(networkPrecCompFake>0,1))+ 1e-4)));

    newG1 = reduce2nNetwork(abs(networkPrecCompFake)>0);
                 
    edgesInNetworkFake(i) = sum(sum(newG1.*triu(origNetwork,1)));
    edgesNotInNetworkFake(i) =  sum(sum(newG1.*triu(~origNetwork,1)));
    toc
end
%% PLOTTING
sensitivity = edgesInNetwork./sum(sum(triu(origNetwork,1)));
falseDiscRate = edgesNotInNetwork./ (edgesInNetwork + edgesNotInNetwork);
figure
subplot(221)
histogram(edgesInNetwork)
title('Edges Found (total = 4)')
subplot(222)
histogram(edgesNotInNetwork)
title('Incorrect Edges Found (total = 6)')
subplot(223)
histogram(corrsNet)
title('Correlation')
subplot(224)
hist([penInComp',penOutComp'])
legend({'Inside','Outside'})
title('Penalization Used')
%%  PLOTTING
figure
subplot(221)
histogram(edgesInNetworkFake)
title('Edges Found (total = 4)')
subplot(222)
histogram(edgesNotInNetworkFake)
title('Incorrect Edges Found (total = 6)')
subplot(223)
histogram(corrsNetFake)
title('Correlation')
subplot(224)
hist([penInCompFake',penOutCompFake'])
legend({'Inside','Outside'})
title('Penalization Used')