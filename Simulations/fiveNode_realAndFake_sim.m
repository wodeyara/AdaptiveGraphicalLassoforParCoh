clear all
S = [6, -.9*exp(1i*pi/4), -.5*exp(1i*pi/8),.55*exp(1i*pi/6), 0;
    -.9*exp(-1i*pi/8), 1.01, 0, 0, -.64*exp(1i*pi/8) ;
    -.5*exp(-1i*pi/8), 0, .6, 0, 0;
    .55*exp(-1i*pi/8), 0, 0, 1, 0;
    0,       -.64*exp(1i*pi/8), 0,0,1.5];
S = .5*(S+S');
S_aug = [S,zeros(5);zeros(5),conj(S)];
realPrec = real2Complex(S_aug,1);

origNetwork = [1,1,1,1,0;
               1,1,0,0,1;
               1,0,1,0,0;
               1,0,0,1,0;
               0,1,0,0,1];
origNetwork = origNetwork - eye(5);

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

GforFit =[(eye(length(origNetwork)) +double(origNetwork)), ... 
    double(origNetwork) ; double(origNetwork), (eye(length(origNetwork)) +double(origNetwork))];
           
lambdasIn = fliplr([.9,.8,.7,.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001]);           
allLambdasOut = fliplr([.9,.8,.7,.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001]); %           
samps = 240;
% SNR = 10;
%%
parfor i = 1:200
    tic
%     noiseAmt = trace(realCov) /(2*length(realCov)* 10 ^ (SNR/10));
%     noiseAmt = noiseAmt * eye(size(realCov));
    data = mvnrnd(zeros(length(realPrec),1),(realPrec)\eye(size(realPrec)),samps)';
    
    reconDataReal = permute(reshape(data', 4,samps/4, size(data,1)),[1,3,2]);

    [networkPrecComp, penInComp(i), penOutComp(i)] = estBestPenalizationQUIC(... 
        reconDataReal, origNetwork,lambdasIn,allLambdasOut, 0);
    
    [networkPrecCompFake, penInCompFake(i), penOutCompFake(i)] = estBestPenalizationQUIC(... 
        reconDataReal, fakeNetwork,lambdasIn,allLambdasOut, 0);    
    
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
%%
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
%% 
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