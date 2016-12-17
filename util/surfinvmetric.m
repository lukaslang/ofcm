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
function [ginv11, ginv12, ginv21, ginv22] = surfinvmetric(Ns, c, rho, xi, detg)
%SURFINVMETRIC Computes the inverse of the metric tensor.
%   
%   [ginv11, ginv12, ginv21, ginv22] = surfinvmetric(Ns, rho, xi, detg) 
%   takes a list of degrees Ns, coefficients c, a radius function rho, 
%   evaluation points xi, and the determinant of the metric, and returns 
%   the metric tensor.
%
%   Ns is a vector of consecutive non-negative integers.
%   c is a vector of length k.
%   rho is a vector of length n.
%   xi is an n-by-2 matrix of coordinates [el, az], where el in [0, pi] and
%   az in [0, 2pi).
%   detg is a vector of length n.
%
%   ginv11, ginv12, ginv21, ginv22 are vectors of length n.

% Get coordinates.
el = xi(:, 1);

% Evaluate first partial derivatives of spherical harmonics.
[d1Yij, d2Yij] = spharmderivn(Ns, xi);

% Compute partial derivatives of rho.
d1rho = d1Yij*c;
d2rho = d2Yij*c;

% Compute metric tensor.
ginv11 = (d2rho.^2 + (rho.^2) .* (sin(el).^2)) ./ detg;
ginv12 = - (d1rho .* d2rho) ./ detg;
ginv21 = ginv12;
ginv22 = (d1rho .^2 + rho.^2) ./ detg;

end