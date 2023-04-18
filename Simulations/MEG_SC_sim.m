% This is code for a simulation to test ability to recover precision
% matrices simulated at the source level, forward modeled to MEG sensors,
% source-localized using a weighted L2 norm inverse before we apply the
% adaptive graphical lasso with the true SC to check ability to recover the
% data. We control the SNR by adding noise at the source level before
% forward modeling the data. 
%%
clear all
% load('C:\Users\wodey\Documents\GitHub\SC_FC_MEG\simulations\dataNeededForSim_MEG.mat');
% get lead field data
load('.\util\dataNeededForSim_MEG.mat')
powerForWeightedL2 = .5;
%following helps set up some vars

run allInverses_MEG_dip.m
samps = 480;

% generate true precision
run genCov_CMVN_SC.m
% % 
% the network for adaptive graphical lasso
GforFit =[(eye(length(origNetwork)) +double(origNetwork)), ... 
    double(origNetwork) ; double(origNetwork), (eye(length(origNetwork)) +double(origNetwork))];
% 
allLambdas = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);           
allLambdasOut = fliplr([.6,.5,.4,.3,.2,.175,.15,.125, .1, .075, .05, .025, .01]);            
%%
samps = 2400; % number of samples to be generated

SNR = 25;

penUsed = 1; % penalization applied to get inverse of lead field

cnter = 1;
% number of iterations:
for ee = 1:200
    run allInverses_MEG_dip.m
    
    % generate new precision
    run genCov_CMVN_SC.m

    tic
    % simulate data
    data =mvnrnd(zeros(length(Q),1),(Q)\eye(size(Q)),samps)';
    
    % store all data:
    allData(ee,:,:) = data;
    allOrigPrec(ee,:,:) = Q;
    
    % generate noise to add to data
    tmp = (Q)\eye(size(Q));
    noiseAmt = trace(tmp) /(length(tmp)* 10 ^ (SNR/10));
    noiseAmtCov = noiseAmt * eye(size(tmp));
    noises = mvnrnd(zeros(length(noiseAmtCov),1),(noiseAmtCov),samps);
    noisesComp = transpose(noises(:,1:114) +1j*noises(:,115:228));
    
    compData = data(1:length(useAreas),:) +1j*data(length(useAreas)+1:2*length(useAreas),:);
    % forward model data to the MEG:
    chanCompData = L*(compData+noisesComp);
    tmp = cat(1, real(chanCompData), imag(chanCompData));
    tmp = cov(tmp');
    
    % define source data:
    sourceData = squeeze(allOpInv(penUsed,:,:))'*(chanCompData);%+noisesComp
    
    sourceDataReal = cat(1,real(sourceData),imag(sourceData));    
    
    sourceDataReal = sourceDataReal*(1/mean(abs(sourceDataReal(:)))); % normalize data
    datareshaped = permute(reshape(sourceDataReal', 4,samps/4, size(sourceDataReal,1)),[1,3,2]);
    datareshaped = datareshaped; %noiseAmt * randn(size(datareshaped));
    
    % apply adaptive graphical lasso
    [networkPrecCompTrue, penInCompTrue(ee,cnter), penOutCompTrue(ee,cnter),~,allDevsReturnTrue(ee,cnter,:,:)] = estBestPenalizationQUIC(... 
        datareshaped, origNetwork,allLambdas,allLambdasOut, 0);
   
    % test recovery of the original precision: 
    corrsOnlySCedges(ee) =corr(Q(triu(GforFit>0,1)), networkPrecCompTrue(triu(GforFit>0,1)));
    
    % this gets the binarized matrix from the augmented binarized precision matrix
    newG1 = reduce2nNetwork(abs(networkPrecCompTrue)>0);
    
    % store all precisions generated
    allNetworksTrue(ee,:,:) = networkPrecCompTrue;  
    
    edgesInNetwork(ee) = sum(sum(newG1.*triu(origNetwork,1)));
    edgesNotInNetwork(ee) =  sum(sum(newG1.*triu(~origNetwork,1)));
    
    sum(sum(newG1.*triu(origNetwork,1))) % how many edges in the true network were recovered?
    sum(sum(newG1.*triu(~origNetwork,1))) % how many edges outside true network?
    cnter = cnter+1;
    toc
end

%% plotting
% sensitivity = edgesIn./sum(sum(triu(fakeNetwork,1)));
% falseDiscRate = edgesOut./ (edgesIn + edgesOut);

% figure
% subplot(221)
% hist([edgesInNetwork'])
% % title('Edges In Fake Network')
% % legend({'In Fake Network', 'Outside True Network'})
% % histogram(sensitivity)
% % title('Sensitivity')
% % title('Sensitiv')
% subplot(222)
% hist([edgesNotInNetwork'])
% % legend({'Outside Fake Network', 'In True Network'})
% % title('Incorrect Edges Found (total = 5721)')
% % histogram(falseDiscRate)
% title('FDR')
% subplot(223)
% histogram(corrsNet)
% title('Correlation')
% subplot(224)
% hist([penInCompTrue',penOutCompTrue'])
% legend({'Inside','Outside'})
% title('Penalization Used')

