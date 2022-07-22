function [network,penalizationIn, penalizationOut,minInd,allDevsReturn] = estBestPenalizationQUIC(data,SC,allLambdas,allLambdasOut, flagForReal)
    % this function takes data (ensembles x sources x samples) and estimates a
    % model using SC-QUIC for each ensemble and estimates (cross validated) deviance on the
    % other ensembles. Using this method it searches for the ideal
    % penalizations (lambda1 and lambda2) to apply to the data, and also the ideal epoch to use to
    % represent the network represented by the data. It does an upper half
    % grid search to identify the ideal penalization (meaning lambda1 >= lambda2).
    % flagForReal is used to indicate whether the data is real/complex
    % valued. If complex valued, a real valued matrix data is expected with
    % the real part occupying the first half of sources and imaginary part
    % the second half of sources. See the following lines of code (dataComplex is complex valued, say from fourier transform
    % and in the ensembles x sources x samples organization): 
    % tmp = real(squeeze(dataComplex(:,:,:)));
    % tmp1 = imag(squeeze(dataComplex(:,:,:)));
    % tmpCompToRealData = cat(2,tmp,tmp1); 
    % For more details refer to Schreir and Scharf, 2010.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% INPUTS
	% data: ensembles x sources x samples matrix
	% SC: the structure of the penalization, named SC for structural connectome 
	% allLambdas: the set of lambda values to pass through, will require some trial and error to determine optimal set, assumed to be sorted
	% flagForReal: 0/1 indicating where the input is real-valued or complex-valued (note this still means the input data 
	%              needs to be real-valued for complex data - real and imaginary being separated)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% OUTPUTS
	% network: the precision estimated from the data for the ensemble with least deviance in CV
	% penalizationIn: lambda1 penalization
	% penalizationOut: lambda2 penalization
	% minInd: ensemble selected
	% allDevsReturn: the deviance values calculated during cross validation
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ani Wodeyar
	% 6/14/2019
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    origNetwork = double(~eye(length(SC))) .* SC > 0;
    if ~flagForReal
%         GforFit =[(eye(length(origNetwork)) +double(origNetwork)),double(origNetwork) ; double(origNetwork), (eye(length(origNetwork)) +double(origNetwork))];
        GforFit =[double(origNetwork),double(origNetwork) ; double(origNetwork), double(origNetwork)];

    else
        GforFit = (eye(length(origNetwork)) +double(origNetwork));
    end
    
    allDevs = zeros(length(allLambdas),length(allLambdasOut),size(data,1),size(data,1));
    for mins = 1:size(data,1)   
        
        parfor lambda = 1:length(allLambdas)
            dataCov = cov(squeeze(data(mins,:,:))');
            scalingVal = max(max(triu(abs(dataCov),1)));
            tmpDev = zeros(length(allLambdasOut),size(data,1));
            cnt = lambda;
            
            for lambdaOut = 1:length(allLambdasOut)
                if allLambdasOut(lambdaOut) < allLambdas(lambda)
                    continue
                end
                offSCedgesLambda = allLambdasOut(lambdaOut);
                [X] = QUIC('default', dataCov, ... 
                    allLambdas(lambda)  *scalingVal *GforFit  +  offSCedgesLambda*scalingVal*  ... 
                    (double(~(GforFit))-eye(length(GforFit))) + min(allLambdas)*scalingVal*eye(length(GforFit)), ...
                     1e-4, 0, 200);
                
                newG = abs(X)>0 ;
                if ~flagForReal
                    newG1 = reduce2nNetwork(newG);
                    GforFit_new = [newG1, newG1.*double(~eye(length(newG1))); newG1.*double(~eye(length(newG1))), newG1];
                else
                    GforFit_new = newG;
                end
                Pweighted = ggmFitHtf(dataCov+ eye(length(dataCov)) *min(allLambdas)  * max(max(triu(abs(dataCov),1))),GforFit_new);
                
                cnt_min = 1;
                dataCovs = zeros(length(setdiff(1:size(data,1),mins)),size(dataCov,2),size(dataCov,2));
                for mins1 = setdiff(1:size(data,1),mins)
                    dataCovs(cnt_min,:,:) = cov(squeeze(data(mins1,:,:))');
                    cnt_min= cnt_min+1;
                end
                for  mins1 = setdiff(1:size(data,1),mins)
                    tmpDev(lambdaOut,mins1) = deviance(squeeze(mean(dataCovs,1)), X);%
                end
                cnt = cnt + 1;
            end
            allDevs(lambda,:,mins,:) = round(tmpDev,2);
        end
    end
    allDevs(allDevs==0) = NaN;
    allDevs(imag(allDevs)~=0) = NaN;
    tmp = (nanmean(squeeze(nanmean(allDevs,4)),3));
    tmp1 = nansum(squeeze(nansum(~isnan(allDevs),4)),3);
	% tries to ensure the penalization chosen has succeeded in creating posDef matrices across all CV
    tmp = tmp.*(tmp1==size(data,1)*(size(data,1)-1)); 
    tmp(tmp==0) = NaN;
    numLocs = size(data,1)*(size(data,1)-1) - 1; 
	% following code ensures we pick some value for the penalization
    while sum(isnan(tmp(:)))== size(tmp,1)^2 && numLocs ~= 0
        tmp = nanmean(squeeze(nanmean(allDevs,4)),3);
        tmp1 = nansum(squeeze(nansum(~isnan(allDevs),4)),3);
        tmp = tmp.*(tmp1>numLocs);
        tmp(tmp==0) = NaN;
        numLocs  = numLocs -1;      
    end

    [~,ind] = min(tmp(:));
    [I1, I2] = ind2sub([size(tmp,1), size(tmp,2)],ind);
    allDevsReturn = (tmp);
    
    [~,minInd] = min(squeeze(nanmean(allDevs(I1,I2,:,:),3)));

    penalizationIn = allLambdas(I1);
    data1 = permute(data,[3,1,2]);
    data1 = reshape(data1,size(data1,1)*size(data1,2),size(data1,3));
    dataCov = cov(data1);
    offSCedgesLambda = allLambdasOut(I2);
    penalizationOut = allLambdasOut(I2);
    X = QUIC('default', dataCov, ... 
        penalizationIn  * max(max(triu(abs(dataCov),1)))*GforFit +  (offSCedgesLambda)*max(max(triu(abs(dataCov),1)))*  ... 
        (double(~(GforFit))-eye(length(GforFit))) + min(allLambdas) * max(max(triu(abs(dataCov),1)))* eye(length(GforFit)), ...
         1e-4, 0, 200);

    newG = abs(X)>0 ;
    if ~flagForReal
        newG1 = reduce2nNetwork(newG);
        GforFit_new = [newG1, newG1.*double(~eye(length(newG1))); newG1.*double(~eye(length(newG1))), newG1];
    else
        GforFit_new = newG;
    end
    
%     network = ggmFitHtf(dataCov+ eye(length(dataCov)) *min(allLambdas)  * max(max(triu(abs(dataCov),1))),GforFit_new);
    network = X;
end

function dev = deviance(S,D)
    %S - covariance matrix of data
    %D - inverse covariance matrix of model
%     logDetD = log(det(S*D));
%     if ~isinf(abs(logDetD))
%         dev = - log(det(D));
%     else
        [~,s,~] = svd(D);
        tmp = (diag(s));
        tmp = tmp(tmp>eps);
        logDetD = sum(log(tmp));
        dev = trace(S*D) - logDetD;
%     end
end
