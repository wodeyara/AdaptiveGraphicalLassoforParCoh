% creating the dipoles for the fsaverage brain
% Algorithm
% Run through all sources, treat the vertices coordinates as vectors and
% take a mean of the vectors to get the centroid. 
% create a matrix containg three orientations for the normal from this
% coordinate in the x, y and z directions (location + orientation is needed
% for writing the dipole file
ROIcentroid=[];
for i = setdiff(unique(labelsTmp),subcortROIs)
    vertsLbls = find(labelsTmp==i);
    
    ROIcentroid(i,:) = mean(Brain.Vertex(vertsLbls, :),1);
end

ROIcentroid = ROIcentroid(useAreas,:);
dipoles_114 = zeros(114*3,6);
for i = 1:114
    dipoles_114(((i-1)*3) +1,:) = [ROIcentroid(i,:), 1,0,0];
    dipoles_114(((i-1)*3) +2,:) = [ROIcentroid(i,:), 0,1,0];
    dipoles_114(((i-1)*3) +3,:) = [ROIcentroid(i,:), 0,0,1];  
end

%% Weighted Centroid
L  =leadFieldCortex;

ROIcentroid=[];
for i = setdiff(unique(labels),subcortROIs)
    vertsLbls = find(labels==i);
    weights = sum(abs(L(grads,vertsLbls).^2),1);
    weights = weights ./ (sum(weights));
    ROIcentroid(i,:) = sum(repmat(weights',1,3) .* Brain.Vertex(vertsLbls, :),1);
end

ROIcentroid = ROIcentroid(useAreas,:);
dipoles_114 = zeros(114*3,6);
for i = 1:114
    dipoles_114(((i-1)*3) +1,:) = [ROIcentroid(i,:), 1,0,0];
    dipoles_114(((i-1)*3) +2,:) = [ROIcentroid(i,:), 0,1,0];
    dipoles_114(((i-1)*3) +3,:) = [ROIcentroid(i,:), 0,0,1];  
end

