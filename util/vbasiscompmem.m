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
function [dy1, dy2] = vbasiscompmem(k, h, x, y, mem)
%VBASISCOMPMEM Evaluates components of vectorial basis functions.
%   
%   [dy1, dy2] = vbasiscomp(k, h, x, y, mem) takes a degree k, a
%   parameter h, a matrix x of locations of the basis functions, a list of 
%   evaluation points y, and returns components of the basis functions 
%   evaluated at y. Function uses up to extra mem bytes for matrix
%   multiplication.
%
%   k is a non-negative integer, h a scalar in (0, 1).
%   x is a m-by-3 matrix where each row is a point on the 2-sphere.
%   y is a n-by-3 matrix where each row is a point on the 2-sphere.
%   mem is a non-negative scalar.
%
%   dy1, dy2 are matrices of size [2*m, n].
m = size(x, 1);
n = size(y, 1);

% Compute junk size.
js = ceil(mem / (8 * m));
assert(js >= 1);

% Compute number of junks.
nj = ceil(n / js);

% Compute dot products between all points.
ip = sparse(m, n);
for p=1:nj
	% Compute indices of junk.
    idx = (p-1)*js+1:min(p*js, n);
    % Compute dot product.
    tmp = x * y(idx, :)';
    % Find indices in support of basis functions.
    [row, col] = find(tmp > h);
    % Compute linear indices.
    tmpidx = sub2ind([m, length(idx)], row, col);
    % Store sparse matrix.
    ip(:, idx) = sparse(row, col, tmp(tmpidx), m, length(idx));
end
[row, col, v] = find(ip);
clear ip;

% Evaluate basis functions at points xi.
b = sparse(row, col, k * ((v - h) .^ (k-1)) / (1 - h)^k, m, n);

% Find coordinates.
[az, el, ~] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
el = pi/2 - el;

% Compute coordinate basis at evaluation points.
[d1, d2] = sphtanbasis([el, az], eye(3));

% Create matrix.
d1b = sparse(m, n);
for p=1:nj
    % Compute indices of junk.
    idx = (p-1)*js+1:min(p*js, n);
    d1b(:, idx) = b(:, idx) .* (x * d1(idx, :)');
end
clear d1;

% Create matrix.
d2b = sparse(m, n);
for p=1:nj
    % Compute indices of junk.
    idx = (p-1)*js+1:min(p*js, n);
    d2b(:, idx) = b(:, idx) .* (x * d2(idx, :)');
end
clear d2;
clear b;

% Compute inverse of metric.
g22 = 1 ./ (sin(el) .^ 2);

% Compute components in the tangent basis.
dy1 = [d1b; bsxfun(@rdivide, d2b, sin(el'))];
dy2 = [bsxfun(@times, g22', d2b); bsxfun(@rdivide, -d1b, sin(el'))];

end