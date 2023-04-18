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
GforFit =[(eye(length(origNetwork)) +double(origNetwork)), ... 
    double(origNetwork) ; double(origNetwork), (eye(length(origNetwork)) +double(origNetwork))];

% penalization parameters: 
lambdasIn = fliplr([.6,.575,.55,.525,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001]);           
allLambdasOut = fliplr([.9,.8,.7,.675,.65,.625, .6,.575,.55,.525,.5,.4,.3,.2,.175,.15]); %           
samps = 24;
SNR = 10;
% test how often we are able to recover the true precision:
parfor i = 1:250
    tic

    % simulate data from precision:
    data = mvnrnd(zeros(length(realPrec),1),(realPrec)\eye(size(realPrec)),samps)';
    
    % separate data into ensembles
    reconDataReal = permute(reshape(data', 4,samps/4, size(data,1)),[1,3,2]);
    
    % estimate precision
    [networkPrecComp, penInComp(i), penOutComp(i)] = estBestPenalizationQUIC(... 
        reconDataReal, origNetwork,lambdasIn,allLambdasOut, 0);
    
    % test recovery of precision:
    % correlation of recovered values:
    corrsNet(i) =(corr((realPrec(triu(GforFit>0,1)|triu(networkPrecComp>0,1)) + 1e-4) ...
                     ,(networkPrecComp(triu(GforFit>0,1)|triu(networkPrecComp>0,1))+ 1e-4)));
    newG1 = reduce2nNetwork(abs(networkPrecComp)>0);
    
    % how many edges:
    edgesInNetwork(i) = sum(sum(newG1.*triu(origNetwork,1)));
    edgesNotInNetwork(i) =  sum(sum(newG1.*triu(~origNetwork,1)));
    toc
end
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
%% so the above works great.
