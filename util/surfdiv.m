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
function v = surfdiv(Ns, cs, xi, bfc, bfcd)
%SURFDIV Evaluates the surface divergence of vectorial basis functions.
%   
%   v = SURFDIV(Ns, cs, xi, bfc, bfcd) takes a surface description Ns, cs,
%   evaluation points xi, components of basis functions bfc, their 
%   derivatives bfcd, and returns the surface divergence of the vector 
%   field spanned by the basis functions.
%
%   Ns is a vector of consecutive non-negative integers.
%   cs is a vector of Fourier coefficients (according to degrees Ns).
%   xi is a [n, 2] matrix of spherical coordinates [el, az] of evaluation
%   points.
%   bfc is a cell array of size [1, 2] each cell of size [2*m, n].
%   bfc is a cell array of size [2, 2] each cell of size [2*m, n].
%
%   v is a matrix of size [2*m, n].

% Compute evaluation points.
Y = sphcoord(xi, eye(3));

% Compute synthesis of evaluation points.
[~, rho] = surfsynth(Ns, Y, cs);

% Evaluate metric.
[g{1,1}, g{1,2}, g{2,1}, g{2,2}] = surfmetric(Ns, cs, rho, xi);

% Compute determinant of metric.
detg = g{1,1} .* g{2,2} - g{1,2} .* g{2,1};

% Compute partial derivatives of metric tensor.
[d1g11, d1g12, d1g21, d1g22, d2g11, d2g12, d2g21, d2g22] = surfmetricderiv(Ns, cs, rho, xi);

d1detg = g{2, 2} .* d1g11 + g{1, 1} .* d1g22 - g{2, 1} .* d1g12 - g{1, 2} .* d1g21;
d2detg = g{2, 2} .* d2g11 + g{1, 1} .* d2g22 - g{2, 1} .* d2g12 - g{1, 2} .* d2g21;

% Compute divergence of each basis function.
v = bfcd{1, 1} + bfcd{2, 2} + bsxfun(@times, bfc{1}, (d1detg ./ (2*detg))') + bsxfun(@times, bfc{2}, (d2detg ./ (2*detg))');

end