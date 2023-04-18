function r = corrOfCorrs(G,H,doPlot)
% given two covariance/precision structures, this function estimates the
% correlation betweek them after removing the bottom triangle of values
% =======================================
% INPUT: 
% G: covariance/precision 1
% H: covariance/precision 2
% doPlot: show a scatterplot of the relationship between G and H
% ========================================
% OUTPUT: 
% r: correlation between G and H
G1 = triu(squeeze(G),1);
G1 = nonzeros(G1(:));

G2 = triu(squeeze(H),1);
G2 = nonzeros(G2(:));

r = corr(G2,G1);

if doPlot ==1
    figure
    scatter(G2,G1)

    title(num2str(r))
end

