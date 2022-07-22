function [inversemat, stat,reconstructed] = inversemodel(gainmat,regmethod,regparam);
%[inversemat, stat] = inversemodel(gainmat,regmethod,regparam);
%Does the inverse model calculation using a L2 norm inverse
%input
%gainmat = gain matrix for a unit source at each location
%regmethod = regularization method - 
%'prctile', sets the percentile of the singular values for the regularization
%'fixed', sets a fixed value of the parameter. 
%
if nargin < 2 
regmethod = 'prctile';
regparam = 50;
end
ndim = length(size(gainmat));
nchan = size(gainmat,1);
nsloc = size(gainmat,2); % number of source locations
if ndim == 3
gainmat = gainmat(:,:);
end;
[u,s,v] = svd(gainmat,'econ');
d = diag(s);
if strcmp(regmethod,'prctile')
regparam = prctile(d,regparam);
end;
d = d./(d.^2+(regparam.^2));
d = diag(d);
inversemat = v*d*u';
if nargout == 3
reconstructed = inversemat*gainmat;
end;
spectrum = inversemat*u;
%spectrum = sum(spectrum.^2,1);
stat.regparam = regparam;
stat.u = u;
stat.s = diag(s);
stat.d = diag(d);
stat.v = v;
stat.spectrum = sum(spectrum.^2,1);


if ndim == 3

inversemat2 = zeros(nsloc,nchan,3);
inversemat2(:,:,1) = inversemat(1:nsloc,:);
inversemat2(:,:,2) = inversemat(nsloc+1:2*nsloc,:);
inversemat2(:,:,3) = inversemat(2*nsloc+1:3*nsloc,:);
inversemat = inversemat2;

temp = zeros(nchan,nsloc,3);
temp(:,:,1) = stat.v(1:nsloc,:)';
temp(:,:,2) = stat.v(nsloc+1:2*nsloc,:)';
temp(:,:,3) = stat.v(2*nsloc+1:3*nsloc,:)';

stat.v = temp;

temp = zeros(3*nsloc,nsloc,3);
temp(:,:,1) = stat.reconstructed(1:nsloc,:)';
temp(:,:,2) = stat.reconstructed(nsloc+1:2*nsloc,:)';
temp(:,:,3) = stat.reconstructed(2*nsloc+1:3*nsloc,:)';

stat.reconstructed = temp;

end;
