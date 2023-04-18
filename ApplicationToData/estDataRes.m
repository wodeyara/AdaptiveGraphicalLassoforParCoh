% Establish weighting structure for lead field: 
W = diag(1./sum((leadField.^2))); %||l_i||^2 summing over channels
W = (W); % based on Hamalainens work from 2006 Neuroimage

% doing a weighted L2 norm inverse: 
opInv = inversemodel(leadField*W,'prctile',10);
sourceData = opInv*double(filtGradData');

% reweighting as this is nice for numerics
sourceData = sourceData*1e12;

% multitaper with 1 taper since we have 1s intervals and I want maximal
% independence between neighboring seconds and also between neighboring
% freq bins
params.tapers = [1,1,1]; % 
params.Fs = 1000;
params.pad = -1;
params.fpass = [1:100];

sourceDataFFT = [];
for i = 1:size(sourceData,1)
    tmp = sourceData(i,:);
    [psds,T, freqs]  = mtspecgramc(tmp,[1,1],params); % don't forget to ensure we get phase from mtspectrumc
    sourceDataFFT(i,:,:) = transpose(psds(:,:));
end

%%
% set up penalization parameters within the connectome:
allLambdas = fliplr([.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001]);           
% set up penalization parameters outside connectome
allLambdasOut = fliplr([.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001]);       

cnt = 1;
% loop over each frequency bin
for i = [1:50]
    tic
    % extract the data for the selected frequency bin
    singleFreqSource = squeeze(sourceDataFFT(:,i,:))';
    % cross-spectrum: 
    data_oneFreq = conj(singleFreqSource) .* singleFreqSource;
    % reduce source data at each source to one dipole using the power at each area:
    cnt1 = 1; dataNew = zeros(size(singleFreqSource,1),size(singleFreqSource,2)/3);
    for j = 3:3:size(leadField,2)
        tmpDipole = data_oneFreq(:,j-2:j);
        [~,~,V] = svd(tmpDipole, 'econ');
        % use the weights on the phase + amp data:
        dataNew(:,cnt1) = singleFreqSource(:,j-2:j)*V(:,1);    
        cnt1 = cnt1 + 1;
    end    
    % reshape data to force 4 ensembles, each with 120 seconds of data
    datareshaped = permute(reshape(dataNew, 4,480/4, size(dataNew,2)),[1,3,2]);
    datareshaped = cat(2,real(datareshaped),imag(datareshaped));
%     datareshaped = datareshaped*(1/mean(abs(datareshaped(:)))); % normalize data
    % estimate precision
    [networkPrecComp, penInComp(cnt), penOutComp(cnt),~,allDevs(cnt,:,:)] = estBestPenalizationQUIC(...
    double(datareshaped), SC>0,allLambdas,allLambdasOut, 0);
    % networkPrecComp is a real-valued matrix represented a complex-valued
    % precision so convert back to complex values: 
    tmp = real2Complex(networkPrecComp,0); % 0 means the matrix was real 
    % remove the part of precision that represents the complementary
    % precision, as we assume propriety for data:
    precComp(cnt,:,:) = abs(tmp(1:114,1:114)).^2./(diag(tmp(1:114,1:114))*diag(tmp(1:114,1:114))');
    tmp1 = squeeze(precComp(cnt,:,:));
    
    % extract data on SC
    precCompSC(cnt,:) = tmp1(find(tril(SC>0,-1)));
    
    % display penalizations selected for data:
    penInComp(cnt)
    penOutComp(cnt)
    cnt = cnt + 1;
    toc
end