% estimate correlation only over edges found in network
% see line 16 for most useful code in file

for k = 1:3
    filesToLoad = {'fakeAndRealNetworks_SNR15_pen10_AGL.mat',...
               'fakeAndRealNetworks_SNR25_pen10_AGL.mat', ...
               'fakeAndRealNetworks_SNR45_pen10_AGL.mat'};
           load(filesToLoad{k})
for i = 1:200
    
     Q = squeeze(allOrigPrec(ee,:,:));
     networkPrecComp = squeeze(allNetworksTrue(i,:,:));             

    % most useful piece of code: this extracts only the correct edges and then
    % estimates correlation between original precision and the estimated
    % precision
     corrsOnlySCedges(i) =corr(Q(triu(GforFit>0,1)), networkPrecComp(triu(GforFit>0,1)));
                 
                 
end
save(filesToLoad{k}, 'corrsOnlySCedges', '-append')
clear
end
%%

for k = 1:4     
for i = 1:200
    
     Q = squeeze(allOrigPrec(i,k,:,:));
     networkPrecComp = squeeze(allNetworksTrue(i,k,:,:));             
    
     corrsOnlySCedges(i,k) =corr(Q(triu(GforFit>0,1)), networkPrecComp(triu(GforFit>0,1)));
       
end
end
save('fakeAndRealNetworks_SNR25_pen35_AGL.mat', 'corrsOnlySCedges', '-append')

%% FOR ONLY SC case
% load('noSL_SC_SNR5to45_sim.mat');

for k = 1:2
for i = 1:200
    
     Q = squeeze(allStartPrec(k,i,:,:));
     networkPrecComp = squeeze(allNetworks(k,i,:,:));             
    
     corrsOnlySCedges(i,k) =corr(Q(triu(GforFit>0,1)), networkPrecComp(triu(GforFit>0,1)));
end
end
