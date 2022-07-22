function z = fisherr2z(r, d, n)
% fisherr2z  - convert Pearson correlation into z values
%
% FORMAT:       z = fisherr2z(r [, inverse [, n]])
%
% Input fields:
%
%       r           correlation value(s) (or z for inverse)
%       inverse     if given and 1x1 logical true, inverse operation
%       n           number of observations in correlation
%
% Output fields:
%
%       z           z values (or r for inverse)
%
% Note: if n given and > 3, z-scores will be (approximately) normally
%       distributed with SD:=1

% Version:  v0.7g
% Build:    9070715
% Date:     Jul-07 2009, 3:55 PM CEST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://wiki.brainvoyager.com/BVQXtools

% argument check
if nargin < 1 || ...
   ~isa(r, 'double') || ...
    any(isinf(r(:)) | isnan(r(:)))
    error( ...
        'BVQXtools:BadArgument', ...
        'Invalid or missing r argument.' ...
    );
end

% correct for SD
if nargin > 2 && ...
    isa(n, 'double') && ...
    numel(n) == 1 && ...
   ~isinf(n) && ...
   ~isnan(n) && ...
    n > 3
    sd = (1 ./ (n - 3)) .^.5;
else
    sd = 1;
end

% compute the desired direction
if nargin < 2 || ...
   ~islogical(d) || ...
   ~d(1)
    z = 0.5 * log((1 + r) ./ (1 - r));
    if ~isreal(z)
        znr = (imag(z) ~= 0);
        z(znr) = Inf * sign(real(z(znr)));
    end
    if sd ~= 1
        z = (1 ./ sd) * z;
    end
else
    if sd ~= 1
        r = sd .* r;
    end
    z = exp(1) .^ (2 .* r);
    z = (z - 1) ./ (z + 1);
end
