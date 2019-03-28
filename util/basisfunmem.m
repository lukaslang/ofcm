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
function b = basisfunmem(k, h, x, y, mem)
%BASISFUNMEM Evaluates basis functions.
%   
%   b = BASISFUNMEM(k, h, x, y, mem) takes a degree k, a parameter h, a
%   matrix x of locations of the basis functions, and a list of evaluation 
%   points y, and returns basis functions evaluated at y. Function uses up
%   to extra mem bytes for matrix multiplication.
%
%   k is a non-negative integer, h a scalar in (0, 1).
%   x is a m-by-3 matrix where each row is a point on the 2-sphere.
%   y is a n-by-3 matrix where each row is a point on the 2-sphere.
%   mem is a non-negative scalar.
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
b = sparse(row, col, ((v - h) .^ k) / (1 - h)^k, m, n);

end