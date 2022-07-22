
i = 1;
namesLeg = {'In Fake Network', 'Outside Fake Network'};
namesLeg2 = {'In True Network', 'Outside True Network'};
colIdx = [ones(200,1),2*ones(200,1)];
g(1,1)=gramm('x',[penInComp(:,i)', penOutComp(:,i)'],'color',namesLeg(colIdx));

%Raw data as scatter plot

g(1,1).stat_bin('edges', [.01,.05:.05:.6] ...
    ); %by default, 'geom' is 'bar', where color groups are side-by-side (dodged)
% g(1,1).set_title('Penalization');
g(1,1).set_names('x','Penalization','y','Count');

g(1,2)=gramm('x',[edgesInFake(:,i)',edgesOutFake(:,i)'],'color',namesLeg(colIdx));

g(1,2).stat_bin('edges', [0:60:720]);
% g(1,2).set_title('');
g(1,2).set_names('x','Edges', 'y','Count');
% g.set_title('');

g(1,3)=gramm('x',[edgesInNetwork(:,i)',edgesOutNetwork(:,i)'],'color',namesLeg2(colIdx));

g(1,3).stat_bin('edges', [0:60:720]);
% g(1,2).set_title('');
g(1,3).set_names('x','Edges', 'y','Count');
% g.set_title('');

g.set_title(['480 samples'])
g.axe_property('Fontsize', 14)

i = 2;
namesLeg = {'In Fake Network', 'Outside Fake Network'};
namesLeg2 = {'In True Network', 'Outside True Network'};
colIdx = [ones(200,1),2*ones(200,1)];
g(2,1)=gramm('x',[penInComp(:,i)', penOutComp(:,i)'],'color',namesLeg(colIdx));

%Raw data as scatter plot

g(2,1).stat_bin('edges', [.01,.05:.05:.6] ...
    ); %by default, 'geom' is 'bar', where color groups are side-by-side (dodged)
% g(1,1).set_title('Penalization');
g(2,1).set_names('x','Penalization','y','Count');

%Boxplots
g(2,2)=gramm('x',[edgesInFake(:,i)',edgesOutFake(:,i)'],'color',namesLeg(colIdx));

g(2,2).stat_bin('edges', [0:60:720]);
% g(1,2).set_title('');
g(2,2).set_names('x','Edges', 'y','Count');
% g.set_title('');

%Boxplots
g(2,3)=gramm('x',[edgesInNetwork(:,i)',edgesOutNetwork(:,i)'],'color',namesLeg2(colIdx));

g(2,3).stat_bin('edges', [0:60:720]);
% g(1,2).set_title('');
g(2,3).set_names('x','Edges', 'y','Count');
% g.set_title('');

g.set_title(['2400 samples'])
g.axe_property('Fontsize', 14)

figure('Position',[100 100 800 550]);

g.draw();

%%

xlims = get(gca,'XLim');
plot(xlims,[720,720],'Linewidth',2, 'Color', [0.64 0.078 0.18])
% xlabel('Samples')
% % ylabel('Edges')
