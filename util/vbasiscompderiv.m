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
function [y1d1, y1d2, y2d1, y2d2] = vbasiscompderiv(k, h, x, y)
%VBASISCOMPDERIV Computes partial derivatives of components of basis
%functions.
%   
%   [y1d1, y1d2, y2d1, y2d2] = vbasiscompderiv(k, h, x, y)
%   takes a degree k, a parameter h, a matrix x of locations of the basis 
%   functions, a list of evaluation points y, and returns components of the
%   basis functions evaluated at y.
%
%   k is a non-negative integer, h a scalar in (0, 1).
%   x is a m-by-3 matrix where each row is a point on the 2-sphere.
%   y is a n-by-3 matrix where each row is a point on the 2-sphere.
%
%   All return values are matrices of size [2*m, n].
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

% Compute second partial derivative of parametrization.
[d11, d12, d21, d22] = sphtanbasisderiv([el, az]);

% Evaluate second derivative of basis functions at points xi.
b2 = sparse(row, col, (k-1) * k * ((v - h) .^ (k-2)) / (1 - h)^k, m, n);

% Compute second partial derivatives of basis functions.
d11b = b2 .* (x * d1') .* (x * d1') + b .* (x * d11');
d21b = b2 .* (x * d2') .* (x * d1') + b .* (x * d21');
d12b = b2 .* (x * d1') .* (x * d2') + b .* (x * d12');
d22b = b2 .* (x * d2') .* (x * d2') + b .* (x * d22');

% Compute partial derivatives of components of basis functions.
y1d1 = [d11b; bsxfun(@rdivide, d12b, sin(el')) - bsxfun(@times, d2b, cos(el') ./ (sin(el').^2))];
y1d2 = [d21b; bsxfun(@rdivide, d22b, sin(el'))];
y2d1 = [-bsxfun(@times, d2b, 2 * cos(el') ./ (sin(el') .^3)) + bsxfun(@rdivide, d12b, sin(el').^2); -bsxfun(@rdivide, d11b, sin(el')) + bsxfun(@times, d1b, cos(el') ./ (sin(el').^2))];
y2d2 = [bsxfun(@rdivide, d22b, sin(el').^2); -bsxfun(@rdivide, d21b, sin(el'))];

end