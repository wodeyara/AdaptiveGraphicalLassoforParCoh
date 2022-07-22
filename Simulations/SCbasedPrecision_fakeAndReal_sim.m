load('C:\Users\wodey\Documents\GitHub\SC_FC_MEG\simulations\dataNeededForSim_MEG.mat');
% clear all

% load('/home/awodeyar/ForGit_MEGpaper/simulations/dataNeededForSim_MEG.mat')
samps = 480;
run genCov_CMVN_SC.m

GforFit =[(eye(length(origNetwork)) +double(origNetwork)), ... 
    double(origNetwork) ; double(origNetwork), (eye(length(origNetwork)) +double(origNetwork))];


allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);           
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);            

SNR = 25; cnter = 1;
for samps = [480,2400]
for ee = 1:200
    tic
    run genCov_CMVN_SC.m
    shuffCol = randperm(length(origNetwork));
    fakeNetwork = origNetwork(shuffCol,shuffCol);

    allOrigNetworks(ee,cnter,:,:) = Q;
    data = mvnrnd(zeros(length(Q),1),(Q)\eye(size(Q)),samps)';
    noiseAmt = trace(covMat) /(length(covMat)* 10 ^ (SNR/10));
    noises = mvnrnd(zeros(length(Q),1),noiseAmt*eye(length(Q)),samps)';
    noises = permute(reshape(noises', 4,samps/4, size(noises,1)),[1,3,2]);
    datareshaped = permute(reshape(data', 4,samps/4, size(data,1)),[1,3,2]);
    datareshaped = datareshaped + noises; %noiseAmt * randn(size(datareshaped));
    
    datareshaped = datareshaped*(1/mean(abs(datareshaped(:)))); % normalize data
        
%     [networkPrecComp, penInCompTrue(ee,cnter), penOutCompTrue(ee,cnter),~,allDevsTrue(ee,cnter,:,:)] = estBestPenalizationQUIC(... 
%         datareshaped, origNetwork,allLambdas,allLambdasOut, 0);
    % if desired to run fake networks do the following:
    [networkPrecComp, penInComp(ee,cnter), penOutComp(ee,cnter),~,allDevs(ee,cnter,:,:)] = estBestPenalizationQUIC(... 
        datareshaped, fakeNetwork,allLambdas,allLambdasOut, 0);
    
    corrsNet(ee,cnter) =(corr((Q(triu(GforFit>0,1)|triu(networkPrecComp~=0,1))) ...
                     ,(networkPrecComp(triu(GforFit>0,1)|triu(networkPrecComp~=0,1)))));
                 
    corrsOnlySCedges(ee,cnter) =corr(Q(triu(GforFit>0,1)), networkPrecComp(triu(GforFit>0,1)));    
    
    allNetworks(ee,cnter,:,:) = networkPrecComp;             
    newG1 = reduce2nNetwork(abs(networkPrecComp)>0);
                 
%     edgesInNetwork(ee,cnter) = sum(sum(newG1.*triu(origNetwork,1)));
%     edgesNotInNetwork(ee,cnter) =  sum(sum(newG1.*triu(~origNetwork,1)));
    edgesIn(ee,cnter) = sum(sum(newG1.*triu(fakeNetwork,1)));
    edgesOut(ee,cnter) = sum(sum(newG1.*triu(double(~fakeNetwork),1)));

    sum(sum(newG1.*triu(fakeNetwork,1)))
     sum(sum(newG1.*triu(double(~fakeNetwork),1)))
    sum(sum(newG1.*triu(origNetwork,1)))
    sum(sum(newG1.*triu(~origNetwork,1)))
    corrsNet(ee)
    
    toc
end
cnter = cnter+1;
end
% sensitivity = edgesInNetwork./sum(sum(triu(origNetwork,1)));
% falseDiscRate = edgesNotInNetwork./ (edgesInNetwork + edgesNotInNetwork);
% 
% figure
% subplot(221)
% histogram(sensitivity)
% title('Sensitivity')
% % title('Edges Found (total = 720)')
% subplot(222)
% % histogram(edgesNotInNetwork)
% % title('Incorrect Edges Found (total = 5721)')
% histogram(falseDiscRate)
% title('FDR')
% subplot(223)
% histogram(corrsNet)
% title('Correlation')
% subplot(224)
% hist([penInComp',penOutComp'])
% legend({'Inside','Outside'})
% title('Penalization Used')
save(['SC_noSL_fake_SNR', num2str(SNR),'480samps.mat'])

%%

% IF desired to run fake network uncomment following: 
% need brain connectivity toolbox (BCT)
% [fakeNetwork] = randomizer_bin_und(origNetwork,1);
% numEdgesToAdd = sum(sum(triu(fakeNetwork.*origNetwork,1))); 
% fakeNetwork = (fakeNetwork - fakeNetwork.*origNetwork); % remove SC edges if any
% tmp1 = double(~fakeNetwork) - eye(length(fakeNetwork)) - origNetwork;
% tmp1 = triu(tmp1,1);
% tmp1Locs = find(tmp1);
% tmp1LocsToUse = randperm(length(tmp1Locs));
% tmp1LocsToUse = tmp1Locs(tmp1LocsToUse(1:numEdgesToAdd));
% fakeNetwork(tmp1LocsToUse) = 1;
% fakeNetwork = double((triu(fakeNetwork,1) + triu(fakeNetwork,1)') > 0);

% shuffCol = randperm(length(origNetwork));
% [V,D] = eig(double(origNetwork));
% fakeNetwork = round(V(shuffCol,:)*D*V(shuffCol,:)');