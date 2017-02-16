% Copyright 2017 Lukas Lang
%
% This file is part of OFCM.
%
%    OFCM is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFCM is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFCM.  If not, see <http://www.gnu.org/licenses/>.
function u = spharmgrad(N, y)
%SPHARMGRAD Computes the gradient of spherical harmonics.
%   
%   u = SPHARMGRAD(N, y) takes a degree N, a list of 
%   evaluation points y, and returns the gradient of spherical harmonics 
%   evaluated at y.
%
%   N is a non-negative integer.
%   y is a n-by-3 matrix where each row is a point on the 2-sphere.
%
%   u is a n-by-m-by-3 matrices, where m are all orders of degree N of
%   spherical harmonics.
assert(N > 0);

% Compute spherical coordinates of evaluation points.
[az, el, ~] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
el = pi/2 - el;

% Compute coordinate basis at evaluation points.
[d1, d2] = sphtanbasis([el, az], eye(3));

% Replicate.
m = 2*N + 1;
d1x = permute(d1, [1, 3, 2]);
d2x = permute(d2, [1, 3, 2]);
d1x = repmat(d1x, 1, m, 1);
d2x = repmat(d2x, 1, m, 1);

% Compute metric.
g22 = repmat(1 ./ (sin(el) .^ 2), 1, m);

% Compute partial derivatives.
[d1, d2] = spharmderiv(N, [el, az]);

% Compute gradient of basis function.
u = bsxfun(@times, d1, d1x) + bsxfun(@times, g22 .* d2, d2x);

end