function [network,networkUnthr] = estImagCoh(data)
% takes data and uses it to define coherence

data = permute(data,[1,3,2]);
data = reshape(data,size(data,1)*size(data,2),size(data,3));

imagcoh = zeros(1000,size(data,2)/2,size(data,2)/2);
parfor i = 1:1000
    dataCov = cov(data(randi(size(data,1),1,size(data,1)),:));
    csd = real2Complex(dataCov,0);
    csd = csd(1:size(data,2)/2,1:size(data,2)/2);
    imagcoh(i,:,:) = abs(imag(csd)./(sqrt(diag(csd)'*diag(csd))));
end

meanImag = squeeze(mean(imagcoh,1));
%fisher transform
imagcoh = .5 * log((1+imagcoh)./(1-imagcoh));
stdImag = reshape(std(reshape(imagcoh,1000,size(imagcoh,2)*size(imagcoh,2)),[],1), ... 
                    size(data,2)/2,size(data,2)/2);
                
lowerBndImag = .5 * log((1+meanImag)./(1-meanImag)) - 1.96*stdImag;
network = meanImag .* (double(lowerBndImag>0) + eye(length(lowerBndImag)));
networkUnthr = meanImag;