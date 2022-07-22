function [network,networkUnthr] = estCoh(data)
% takes data and uses it to define coherence


data = permute(data,[1,3,2]);
data = reshape(data,size(data,1)*size(data,2),size(data,3));

coh = zeros(1000,size(data,2)/2,size(data,2)/2);
parfor i = 1:1000
    dataCov = cov(data(randi(size(data,1),1,size(data,1)),:));
    csd = real2Complex(dataCov,0);
    csd = csd(1:size(data,2)/2,1:size(data,2)/2);
    coh(i,:,:) = abs((csd)./(sqrt(diag(csd)'*diag(csd))));
end

meanCoh = squeeze(mean(coh,1));
%fisher transform
coh = .5 * log((1+coh)./(1-coh));
stdImag = reshape(std(reshape(coh,1000,size(coh,2)*size(coh,2)),[],1), ... 
                    size(data,2)/2,size(data,2)/2);
                
lowerBndImag = .5 * log((1+meanCoh)./(1-meanCoh)) - 1.96*stdImag;
thr = triu(double(lowerBndImag>0),1) ;
network = meanCoh .* (thr + thr' + eye(length(lowerBndImag)));
networkUnthr = meanCoh;