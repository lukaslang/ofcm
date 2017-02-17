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
function [y1d1, y1d2, y2d1, y2d2] = vbasiscompderivmem(k, h, x, y, mem)
%VBASISCOMPDERIVMEM Computes partial derivatives of components of basis
%functions.
%   
%   [y1d1, y1d2, y2d1, y2d2] = vbasiscompderiv(k, h, x, y)
%   takes a degree k, a parameter h, a matrix x of locations of the basis 
%   functions, a list of evaluation points y, and returns components of the
%   basis functions evaluated at y. Uses up to mem bytes for matrix
%   multiplication.
%
%   k is a non-negative integer, h a scalar in (0, 1).
%   x is a m-by-3 matrix where each row is a point on the 2-sphere.
%   y is a n-by-3 matrix where each row is a point on the 2-sphere.
%
%   All return values are matrices of size [2*m, n].
m = size(x, 1);
n = size(y, 1);

% Compute junk size.
js = ceil(mem / (3 * 8 * m));
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

% Find coordinates.
[az, el, ~] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
el = pi/2 - el;

% Compute coordinate basis at evaluation points.
[d1, d2] = sphtanbasis([el, az], eye(3));

% Compute second partial derivative of parametrization.
[d11, d12, d21, d22] = sphtanbasisderiv([el, az]);

% Evaluate basis functions at points xi.
b = sparse(row, col, k * ((v - h) .^ (k-1)) / (1 - h)^k, m, n);

% Evaluate second derivative of basis functions at points xi.
b2 = sparse(row, col, (k-1) * k * ((v - h) .^ (k-2)) / (1 - h)^k, m, n);

% Create matrix.
d11b = sparse(m, n);
for p=1:nj
    % Compute indices of junk.
    idx = (p-1)*js+1:min(p*js, n);
    d11b(:, idx) = b2(:, idx) .* (x * d1(idx, :)') .* (x * d1(idx, :)') + b(:, idx) .* (x * d11(idx, :)');
end

% Create matrix.
d12b = sparse(m, n);
for p=1:nj
    % Compute indices of junk.
    idx = (p-1)*js+1:min(p*js, n);
    d12b(:, idx) = b2(:, idx) .* (x * d1(idx, :)') .* (x * d2(idx, :)') + b(:, idx) .* (x * d12(idx, :)');
end

% Create matrix.
d21b = sparse(m, n);
for p=1:nj
    % Compute indices of junk.
    idx = (p-1)*js+1:min(p*js, n);
    d21b(:, idx) = b2(:, idx) .* (x * d2(idx, :)') .* (x * d1(idx, :)') + b(:, idx) .* (x * d21(idx, :)');
end

% Create matrix.
d21b = sparse(m, n);
for p=1:nj
    % Compute indices of junk.
    idx = (p-1)*js+1:min(p*js, n);
    d22b(:, idx) = b2(:, idx) .* (x * d2(idx, :)') .* (x * d2(idx, :)') + b(:, idx) .* (x * d22(idx, :)');
end
clear d11;
clear d12;
clear d21;
clear d22;

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

% Compute partial derivatives of components of basis functions.
y1d1 = [d11b; bsxfun(@rdivide, d12b, sin(el')) - bsxfun(@times, d2b, cos(el') ./ (sin(el').^2))];
y1d2 = [d21b; bsxfun(@rdivide, d22b, sin(el'))];
y2d1 = [-bsxfun(@times, d2b, 2 * cos(el') ./ (sin(el') .^3)) + bsxfun(@rdivide, d12b, sin(el').^2); -bsxfun(@rdivide, d11b, sin(el')) + bsxfun(@times, d1b, cos(el') ./ (sin(el').^2))];
y2d2 = [bsxfun(@rdivide, d22b, sin(el').^2); -bsxfun(@rdivide, d21b, sin(el'))];

end