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
function [d1g11, d1g12, d1g21, d1g22, d2g11, d2g12, d2g21, d2g22] = surfmetricderiv(Ns, c, rho, xi)
%SURFMETRICDERIV Computes partial derivatives of the metric tensor.
%   
%   [d1g11, d1g12, d1g21, d1g22, d2g11, d2g12, d2g21, d2g22] = surfmetricderiv(Ns, rho, xi)
%   takes list of degrees Ns, coefficients c, a radius function rho and 
%   evaluation points xi, and returns partial derivatives of the 
%   metric tensor of a sphere-like surface.
%
%   Ns is a vector of consecutive non-negative integers.
%   c is a vector of length k.
%   rho is a vector of length n.
%   xi is an n-by-2 matrix of coordinates [el, az], where el in [0, pi] and
%   az in [0, 2pi).
%
%   Output arguments are vectors of length n.

% Get coordinates.
el = xi(:, 1);

% Evaluate first partial derivatives of spherical harmonics.
[d1Yij, d2Yij] = spharmderivn(Ns, xi);

% Evaluate second partial derivatives of spherical harmonics.
[d11Yij, d12Yij, d21Yij, d22Yij] = spharmderivn2(Ns, xi);

% Compute partial derivatives of rho.
d1rho = d1Yij*c;
d2rho = d2Yij*c;

% Compute second partial derivatives of rho.
d11rho = d11Yij*c;
d12rho = d12Yij*c;
d21rho = d21Yij*c;
d22rho = d22Yij*c;

% Compute partial derivatives with respect to elevation.
d1g11 = 2 * d11rho .* d1rho + 2 * d1rho;
d1g12 = d11rho .* d2rho + d1rho .* d12rho;
d1g21 = d1g12;
d1g22 = 2 * d12rho .* d2rho + 2 * d1rho .* (sin(el).^2) + (rho .^2) .* sin(2*el);

% Compute partial derivatives with respect to azimuthal.
d2g11 = 2 * d21rho .* d1rho + 2 * d2rho;
d2g12 = d21rho .* d2rho + d1rho .* d22rho;
d2g21 = d2g12;
d2g22 = 2 * d22rho .* d2rho + 2 * d2rho .* (sin(el).^2);

end