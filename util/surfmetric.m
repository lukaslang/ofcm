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
function [g11, g12, g21, g22] = surfmetric(Ns, c, rho, xi)
%SURFMETRIC Computes the metric tensor of a sphere-like surface.
%   
%   [g11, g12, g21, g22] = surfmetric(Ns, rho, xi) takes list of degrees 
%   Ns, coefficients c, a radius function rho and evaluation points xi, and
%   returns the metric tensor.
%
%   Ns is a vector of consecutive non-negative integers.
%   c is a vector of length k.
%   rho is a vector of length n.
%   xi is an n-by-2 matrix of coordinates [el, az], where el in [0, pi] and
%   az in [0, 2pi).
%
%   g11, g12, g21, g22 are vectors of length n.

% Get coordinates.
el = xi(:, 1);

% Evaluate first partial derivatives of spherical harmonics.
[d1Yij, d2Yij] = spharmderivn(Ns, xi);

% Compute partial derivatives of rho.
d1rho = d1Yij*c;
d2rho = d2Yij*c;

% Compute metric tensor.
g11 = d1rho.^2 + rho.^2;
g12 = d1rho .* d2rho;
g21 = g12;
g22 = d2rho .^2 + (rho.^2) .* (sin(el).^2);

end