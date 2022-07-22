%% Inverse solution
% this script takes the already loaded leadfield - should be called
% red_L_SVD and the labels for areas for all vertices and returns inverses
% for all percentiles from 1:91 in increments of 10 in the matrix allInv. 

numSources = size(leadField,2)/3;

% generate random weighting matrix for lead field every run
theta = rand(114,1) * (2*pi-1e-6) + 1e-6; %[1e-6,2*pi]
u = rand(114,1)*2 - 1; %[-1,1]

Xfull = [u,theta]';
clear XfullCart;
[XfullCart(:,1),XfullCart(:,2),XfullCart(:,3)] =  ...
    sph2cart(Xfull(2,:)', acos(Xfull(1,:))', ones(numSources,1));
XfullCart = reshape(XfullCart',342,1);
weightingMatrix = zeros(size(leadField,2),numSources);
for i = 3:3:342;
   weightingMatrix((1:3)+((i/3)-1)*3,i/3) = XfullCart((i-2) :i);
end

L = leadField * weightingMatrix; % channels x (dipoles*3) * weighting matrix is (dipoles*3) x channels

% Currently only the weighted L2:
goodchans = 1:size(L,1);
W = 1./diag(sum((L.^2))); %||l_i||^2 summing over channels
W(isinf(W)) = 0;
W = (W.^powerForWeightedL2);
% the following step reweights the Weights so that they don't change the
% diagonal of lead field to 1s
% W = W * (1/mean(nonzeros(W)));

cnt =1;
allSourcePens = [1:91];
allOpInv = zeros(length(allSourcePens), size(L,1), size(L,2));
for penalization=   allSourcePens  
    [opInv, stat,reconstructed] = inversemodel(L*W,'prctile',penalization );
    opInv = opInv';
    allOpInv(cnt , : , :) = opInv;
    cnt = cnt +1;
end

clear W opInv cnt weightingMatrix theta u XfullCart