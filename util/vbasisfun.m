% Copyright 2015 Lukas Lang
%
% This file is part of OFDM.
%
%    OFDM is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFDM is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFDM.  If not, see <http://www.gnu.org/licenses/>.
function u = vbasisfun(k, h, x, y)
%VBASISFUN Evaluates vectorial basis functions.
%   
%   [u, v] = VBASISFUN(k, h, x, y) takes a degree k, a parameter h, a 
%   matrix x of locations of the basis functions, a list of evaluation 
%   points y, and returns basis functions evaluated at y.
%
%   k is a non-negative integer, h a scalar in (0, 1).
%   x is a m-by-3 matrix where each row is a point on the 2-sphere.
%   y is a n-by-3 matrix where each row is a point on the 2-sphere.
%
%   u is of size [2*m, n, 3] matrices.

% Compute components of basis functions.
[dy1, dy2] = vbasiscomp(k, h, x, y);

% Find coordinates.
[az, el, ~] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
el = pi/2 - el;

% Compute coordinate basis at evaluation points.
[d1, d2] = sphtanbasis([el, az], eye(3));

% Compute gradient of basis function and orthogonal to gradient.
u = bsxfun(@times, full(dy1), permute(d1, [3, 1, 2])) + bsxfun(@times, full(dy2), permute(d2, [3, 1, 2]));

end