startVec = [0 0 0 0 .25 .25 .25 .25 zeros(1,808)];
for i = 1:102 
    transformMatMag(i,:) = startVec;
    startVec = circshift(startVec',8)';    
end

magsLeadField = transformMatMag*linop;

leadFieldTmp = zeros(510,size(linop,2));
leadFieldTmp(5:5:510,:) = magsLeadField;
leadFieldTmp(setdiff(1:510,5:5:510), :) = linop(~sum(transformMatMag,1),:);

leadField = gradTransformMatrix*leadFieldTmp;