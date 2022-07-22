map = (brewermap(14,'Spectral'));
i = 1;
load('SC_SLpen1_SNR25_480samps_2400_compAlgs.mat', 'edgesInNetwork')
load('SC_SLpen1_SNR25_480samps_2400_compAlgs.mat', 'edgesNotInNetwork')
newEdges(:,1:2:5) = squeeze(edgesInNetwork(:,i,:));
newEdges(:,2:2:6) = squeeze(edgesNotInNetwork(:,i,:));
load('SC_SLpen1_SNR25_480samps_2400_compAlgs.mat', 'corrsOnlySCedges')
corrsSC(:,1:3) = squeeze(corrsOnlySCedges(:,i,:));

load('SNR25_noiseAtSource_Samps480_2400.mat', 'edgesInNetwork')
load('SNR25_noiseAtSource_Samps480_2400.mat', 'edgesNotInNetwork')
newEdges(:,7) = squeeze(edgesInNetwork(:,i));
newEdges(:,8) = squeeze(edgesNotInNetwork(:,i));
load('SNR25_noiseAtSource_Samps480_2400.mat', 'corrsOnlySCedges')
corrsSC(:,4) = squeeze(corrsOnlySCedges(:,i));

figure(1)
subplot(221);
violinplot(newEdges,nonzeros((diag([1:8])*ones(8,200))'),map([3,4,6,7,9,10,12,13],:))
set(gca,'YScale','log')
set(gca,'XTick', [1.5,3.5,5.5,7.5])
set(gca,'Fontsize', 14)
set(gca,'XTickLabel',{'L2 norm', 'Coherence', 'Imag Coh','AGL'})
ylabel('Edges')
view([90,90])
grid on
hold on
plot([0,10], [720,720],'Linewidth',2, 'Color', [0.64 0.078 0.18])

subplot(222)
h = violinplot(corrsSC,[ones(1,200),2*ones(1,200),3*ones(1,200),4*ones(1,200)],map(3:3:end,:))
set(gca,'Fontsize', 14)
set(gca,'XTickLabel',{'L2 norm', 'Coherence', 'Imag Coh','AGL'})
ylabel('Correlation')
% legend({'L2 norm', 'Coherence', 'Imag Coh'})
view([90,90])
grid on
