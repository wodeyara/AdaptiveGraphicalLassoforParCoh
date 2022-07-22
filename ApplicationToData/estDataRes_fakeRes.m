%%

allLambdas = fliplr([.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001]);           
allLambdasOut = fliplr([.2,.175,.15,.125, .1, .075, .05, .025, .01,.005,.001]);       

origNetwork = SC > 0;
shuffCol = randperm(length(origNetwork));
fakeNetwork = origNetwork(shuffCol,shuffCol);

cnt = 1;
for i = [1:50]
    tic
    singleFreqSource = squeeze(sourceDataFFT(:,i,:))';
    data_oneFreq = conj(singleFreqSource) .* singleFreqSource;
    % reduce xyz to one dipole using the power at each area:
    cnt1 = 1; dataNew = zeros(size(singleFreqSource,1),size(singleFreqSource,2)/3);
    for j = 3:3:size(leadField,2)
        tmpDipole = data_oneFreq(:,j-2:j);
        [~,~,V] = svd(tmpDipole, 'econ');
        % use the weights on the phase + amp data:
        dataNew(:,cnt1) = singleFreqSource(:,j-2:j)*V(:,1);    
        cnt1 = cnt1 + 1;
    end    
    
    datareshaped = permute(reshape(dataNew, 4,480/4, size(dataNew,2)),[1,3,2]);
    datareshaped = cat(2,real(datareshaped),imag(datareshaped));
%     datareshaped = datareshaped*(1/mean(abs(datareshaped(:)))); % normalize data

    [networkPrecComp, penInComp(cnt), penOutComp(cnt),~,allDevs(cnt,:,:)] = estBestPenalizationQUIC(...
    double(datareshaped), fakeNetwork,allLambdas,allLambdasOut, 0);
    tmp = real2Complex(networkPrecComp,0);
    precComp(cnt,:,:) = abs(tmp(1:114,1:114)).^2./(diag(tmp(1:114,1:114))*diag(tmp(1:114,1:114))');
    tmp1 = squeeze(precComp(cnt,:,:));
    precCompSC(cnt,:) = tmp1(find(tril(SC>0,-1)));
    penInComp(cnt)
    penOutComp(cnt)
    cnt = cnt + 1;
    toc
end