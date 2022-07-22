function r = corrOfCorrs(G,H,doPlot)

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

