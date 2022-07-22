clear all
load('.\util\dataNeededForSim_MEG.mat')

% load('/home/awodeyar/ForGit_MEGpaper/simulations/dataNeededForSim_MEG.mat')
SNR = 25;
samps = 480;
run genCov_CMVN_SC.m

GforFit =[(eye(length(origNetwork)) +double(origNetwork)), ... 
    double(origNetwork) ; double(origNetwork), (eye(length(origNetwork)) +double(origNetwork))];

allLambdasOut = fliplr([.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001,.0005,.0001,.00005,.00001]);    

SNR = 25;cnter =1;
for samps = [480, 2400]
for ee = 1:200
    tic
    run genCov_CMVN_SC.m
    
    allOrigNetworks(ee,cnter,:,:) = Q;
    data = mvnrnd(zeros(length(Q),1),(Q)\eye(size(Q)),samps)';
    noiseAmt = trace(covMat) /(length(covMat)* 10 ^ (SNR/10));
    noises = mvnrnd(zeros(length(Q),1),noiseAmt*eye(length(Q)),samps)';
    noises = permute(reshape(noises', 4,samps/4, size(noises,1)),[1,3,2]);
    datareshaped = permute(reshape(data', 4,samps/4, size(data,1)),[1,3,2]);
    datareshaped = datareshaped + noises; %noiseAmt * randn(size(datareshaped));
   
    
    datareshaped = datareshaped*(1/mean(abs(datareshaped(:)))); % normalize data
        
    [netTmp,netTmp2,penalizationIn(ee,cnter), thresholdOut(ee,cnter),minInd,allDevsReturn(ee,cnter,:,:)] = estBestPenalization_L2...
                (datareshaped, SC(useAreas,useAreas)>0, allLambdasOut, 0); 
    [allNetworks(ee,cnter,2,:,:),allNetworksUnthr(ee,cnter,2,:,:)] = estCoh(datareshaped);
    [allNetworks(ee,cnter,3,:,:),allNetworksUnthr(ee,cnter,3,:,:)] = estImagCoh(datareshaped);
    
    netTmp = real2Complex(netTmp,0);
    netTmp = abs(netTmp(1:114,1:114))/sqrt(diag(netTmp(1:114,1:114))'*diag(netTmp(1:114,1:114)));
    allNetworks(ee,cnter,1,:,:) = netTmp;
    
    netTmp2 = real2Complex(netTmp2,0);
    netTmp2 = abs(netTmp2(1:114,1:114))/sqrt(diag(netTmp2(1:114,1:114))'*diag(netTmp2(1:114,1:114)));
    allNetworksUnthr(ee,cnter,1,:,:) = netTmp2;
    
    for kk = 1:3
        net_red = squeeze(allNetworks(ee,cnter,kk,:,:));
       
        newG1 = abs(net_red)>0;
        edgesInNetwork(ee,cnter,kk) = sum(sum(newG1.*triu(origNetwork,1)));
        edgesNotInNetwork(ee,cnter,kk) =  sum(sum(newG1.*triu(~origNetwork,1)));

        sum(sum(newG1.*triu(origNetwork,1)))
        sum(sum(newG1.*triu(~origNetwork,1)))
    end
    toc
end
cnter=cnter+1;
end
% sensitivity = edgesInNetwork./sum(sum(triu(origNetwork,1)));
% falseDiscRate = edgesNotInNetwork./ (edgesInNetwork + edgesNotInNetwork);
% 

% save(['SC_noSL_SNR', num2str(SNR),'.mat'])
%%
% compute correlation correctly
for ee = 1:200
    for cnter = 1:2
    Q1 = real2Complex(squeeze(allOrigNetworks(ee,cnter,:,:)),0);
	Q1 = Q1(1:114,1:114);
    Q1 = abs(Q1)./(sqrt(diag(Q1)'*diag(Q1)));
    
    for kk = 1:3
        net_red = squeeze(allNetworks(ee,cnter,kk,:,:));
        corrsOnlySCedges(ee,cnter,kk) = corr(Q1(triu(SC(useAreas,useAreas)>0,1)),...
            net_red(triu(SC(useAreas,useAreas)>0,1)));

    end
    end
end

save(['SC_noSL_SNR', num2str(SNR),'.mat'])