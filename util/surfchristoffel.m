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
function [c111, c112, c121, c122, c221, c222] = surfchristoffel(Ns, c, rho, xi)
%SURFCHRISTOFFEL Computes the Christoffel symbols of a sphere-like surface.
%   
%   [c111, c112, c121, c122, c221, c222] = surfchristoffel(Ns, c, rho, xi)
%   takes list of degrees Ns, coefficients c, a radius function rho and 
%   evaluation points xi, and returns the Christoffel symbols.
%
%   Ns is a vector of consecutive non-negative integers.
%   c is a vector of length k.
%   rho is a vector of length n.
%   xi is an n-by-2 matrix of coordinates [el, az], where el in [0, pi] and
%   az in [0, 2pi).
%
%   All output arguments are vectors of length n.

% Evaluate metric.
[g11, g12, g21, g22] = surfmetric(Ns, c, rho, xi);

% Compute determinant of metric.
detg = g11 .* g22 - g12 .* g21;

% Evaluate inverse of metric.
[ginv11, ginv12, ginv21, ginv22] = surfinvmetric(Ns, c, rho, xi, detg);

% Compute partial derivatives of metric tensor.
[d1g11, ~, d1g21, d1g22, d2g11, d2g12, ~, d2g22] = surfmetricderiv(Ns, c, rho, xi);

% Compute Christoffel symbols.
c111 = ginv12 .* d1g21 + 0.5 * (ginv11 .* d1g11 - ginv12 .* d2g11);
c112 = ginv22 .* d1g21 + 0.5 * (ginv21 .* d1g11 - ginv22 .* d2g11);
c121 = 0.5 * (ginv11 .* d2g11 + ginv12 .* d1g22);
c122 = 0.5 * (ginv21 .* d2g11 + ginv22 .* d1g22);
c221 = ginv11 .* d2g12 + 0.5 * (ginv12 .* d2g22 - ginv11 .* d1g22);
c222 = ginv21 .* d2g12 + 0.5 * (ginv22 .* d2g22 - ginv21 .* d1g22);

end