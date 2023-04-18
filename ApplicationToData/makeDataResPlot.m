%% This code plots the results from estDataRes

% We will use only those frequencies when the SC was considered useful.
% First how many frequencies showed this (as a ratio of the frequency bands
% they appear in):
% Result 1:
tmpDiff = penOutComp - penInComp;
sum(tmpDiff(1:3) > 0)/3
sum(tmpDiff(4:7) > 0)/4
sum(tmpDiff(8:13) > 0)/6
sum(tmpDiff(14:29) > 0)/15
sum(tmpDiff(30:50) > 0)/21

% Result 2: Average across the times when it was true that the SC was
% useful and show me:

tmpDiffInds = find(tmpDiff>0);
redPrecCoh = precComp(tmpDiffInds,:,:);
[~,tmpMaxInd] = max(redPrecCoh,[],1);
tmpMaxInd = squeeze(tmpMaxInd).*(SC>0) .* double(~eye(114));
tmpMaxInd(tmpMaxInd>0) = tmpDiffInds(tmpMaxInd(tmpMaxInd>0));
indsSC = find(tril(SC,-1));

% find the maximum edges in each and plot those, and only the ones on the
% SC:
subplot(131)
deltaPrecCoh = zeros(114,114);
for i = 1:length(indsSC)
    t = tmpMaxInd(indsSC(i));
    if t>0 && t<4
        t1 = squeeze(precComp(t,:,:));
        deltaPrecCoh(indsSC(i)) = t1(indsSC(i));
    end
end
imagesc((deltaPrecCoh + deltaPrecCoh')>0)
caxis([0,.15])
% colormap((brewermap([],'PuOr')))
axis square

subplot(132)
thetaPrecCoh = zeros(114,114);
for i = 1:length(indsSC)
    t = tmpMaxInd(indsSC(i));
    if t>3 && t<8
        t1 = squeeze(precComp(t,:,:));
        thetaPrecCoh(indsSC(i)) = t1(indsSC(i));
    end
end
imagesc((thetaPrecCoh + thetaPrecCoh')>0)
caxis([0,.15])
% colormap((brewermap([],'PuOr')))
axis square

subplot(133)
betaPrecCoh = zeros(114,114);
for i = 1:length(indsSC)
    t = tmpMaxInd(indsSC(i));
    if t>13 && t<30
        t1 = squeeze(precComp(t,:,:));
        betaPrecCoh(indsSC(i)) = t1(indsSC(i));
    end
end
imagesc((betaPrecCoh + betaPrecCoh')>0)
caxis([0,.15])
% colormap((brewermap([],'PuOr')))
axis square

subplot(131); title('Delta')
set(gca,'Fontsize', 14)
subplot(132); title('Theta')
set(gca,'Fontsize', 14)
subplot(133); title('Beta')
set(gca,'Fontsize', 14)

%%

% imagesc((thetaPrecCoh + thetaPrecCoh')>0)
% caxis([0,.15])
% % colormap((brewermap([],'PuOr')))
% axis square
% set(gca,'XTick', 1:114, 'XTickLabel', roiNames_60,'Ytick', 1:114, 'YTickLabel', roiNames_60)
% xtickangle(45)
