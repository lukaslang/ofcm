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
function [dy1, dy2] = vbasiscomp(k, h, x, y)
%VBASISCOMP Evaluates components of vectorial basis functions.
%   
%   [dy1, dy2] = vbasiscomp(k, h, x, y) takes a degree k, a
%   parameter h, a matrix x of locations of the basis functions, a list of 
%   evaluation points y, and returns components of the basis functions 
%   evaluated at y.
%
%   k is a non-negative integer, h a scalar in (0, 1).
%   x is a m-by-3 matrix where each row is a point on the 2-sphere.
%   y is a n-by-3 matrix where each row is a point on the 2-sphere.
%
%   dy1, dy2 are vectors of length 2*n.
m = size(x, 1);
n = size(y, 1);

% Compute dot products between all points.
ip = x*(y');
idx = find(ip > h);
v = ip(idx);
[row, col] = ind2sub([m, n], idx);

% Evaluate basis functions at points xi.
b = sparse(row, col, k * ((v - h) .^ (k-1)) / (1 - h)^k, m, n);

% Find coordinates.
[az, el, ~] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
el = pi/2 - el;

% Compute coordinate basis at evaluation points.
[d1, d2] = sphtanbasis([el, az], eye(3));

% Compute partial derivatives of basis function.
d1b = b .* (x * d1');
d2b = b .* (x * d2');
clear b;
clear d1;
clear d2;

% Compute inverse of metric.
g22 = 1 ./ (sin(el) .^ 2);

% Compute components in the tangent basis.
dy1 = [d1b; bsxfun(@rdivide, d2b, sin(el'))];
dy2 = [bsxfun(@times, g22', d2b); bsxfun(@rdivide, -d1b, sin(el'))];

end