for i = [4]

clear g

namesLeg = {'In Network', 'Outside Network'};
colIdx = [ones(200,1),2*ones(200,1)];
g(1,1)=gramm('x',[penInComp(i,:)', penOutComp(i,:)'],'color',namesLeg(colIdx));

%Raw data as scatter plot

g(1,1).stat_bin('edges', [.01,.05:.05:.6] ...
    ); %by default, 'geom' is 'bar', where color groups are side-by-side (dodged)
% g(1,1).set_title('Penalization');
g(1,1).set_names('x','Penalization','y','Count');

%Boxplots
g(1,2)=gramm('x',[edgesInNetwork(i,:)',edgesNotInNetwork(i,:)'],'color',namesLeg(colIdx));

g(1,2).stat_bin('edges', [0:60:720]);
% g(1,2).set_title('');
g(1,2).set_names('x','Edges', 'y','Count');
% g.set_title('');

g(1,3) = gramm('x', corrsOnlySCedges(:,i), 'color', ones(200,1));
g(1,3).stat_bin('edges',[0:.05:1])
g(1,3).set_names('x', 'Correlation', 'y','Count');
g.set_title(['SNR = ', num2str(25)])
g.axe_property('Fontsize', 14)
figure('Position',[100 100 800 550]);

g.draw();
end


%%


clear g

namesLeg = {'In Network', 'Outside Network'};
colIdx = [ones(200,1),2*ones(200,1)];
g(1,1)=gramm('x',[penInComp', penOutComp'],'color',namesLeg(colIdx));

%Raw data as scatter plot

g(1,1).stat_bin('edges', [.01,.05:.05:.6] ...
    ); %by default, 'geom' is 'bar', where color groups are side-by-side (dodged)
% g(1,1).set_title('Penalization');
g(1,1).set_names('x','Penalization','y','Count');

%Boxplots
g(1,2)=gramm('x',[edgesInNetwork',edgesNotInNetwork'],'color',namesLeg(colIdx));

g(1,2).stat_bin('edges', [0:10:150]);
% g(1,2).set_title('');
g(1,2).set_names('x','Edges', 'y','Count');
% g.set_title('');

g.set_title(['SNR = ', num2str(SNR)])
g.axe_property('Fontsize', 14)
figure('Position',[100 100 800 550]);

g.draw();