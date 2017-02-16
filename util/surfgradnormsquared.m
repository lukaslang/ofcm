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
function f = surfgradnormsquared(Ns, c, xi)
%SURFGRADNORMSQUARED Computes the squared norm of the surface gradient of a
%function on the sphere given in spherical harmonics expansion.
%   
%   f = surfgradnormsquared(Ns, c, xi) takes list of degrees Ns, 
%   coefficients c, a radius function rho and evaluation points xi, and
%   returns the squared norm of the surface gradient.
%
%   Ns is a vector of consecutive non-negative integers.
%   c is a vector of length k.
%   rho is a vector of length n.
%   xi is an n-by-2 matrix of coordinates [el, az], where el in [0, pi] and
%   az in [0, 2pi).
%
%   f is a vector of length n.

% Get coordinates.
el = xi(:, 1);

% Evaluate first partial derivatives of spherical harmonics.
[d1Yij, d2Yij] = spharmderivn(Ns, xi);

% Compute partial derivatives of rho.
d1rho = d1Yij*c;
d2rho = d2Yij*c;

f = d1rho.^2 + (d2rho .^2) ./ (sin(el).^2);

end