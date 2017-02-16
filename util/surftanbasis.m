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
function [d1, d2] = surftanbasis(Ns, cs, xi)
%SURFTANBASIS Computes tangent basis on surface.
%   
%   [d1, d2] = SURFTANBASIS(cs, xi) takes degrees Ns and corresponding 
%   Fourier coefficients cs and evaluation points xi and returns the
%   tangent basis.
%
%   Ns is a list os non-negative consecutive integers.
%   cs is a vector of Fourier coefficients.
%   xi is an n-by-2 matrix of coordinates [el, az], where el in [0, pi] and
%   az in [0, 2pi).
%
%   d1 and d2 are n-by-3 matrices.

% Get coordinates.
el = xi(:, 1);
az = xi(:, 2);

% Compute tangent basis of sphere.
[d1s, d2s] = sphtanbasis(xi, eye(3));

% Compute points on sphere.
x = [sin(el).*cos(az), sin(el).*sin(az), cos(el)];

% Compute synthesis of evaluation points.
[~, rho] = surfsynth(Ns, x, cs);

% Compute partial derivatives of spherical harmonics.
[d1Y, d2Y] = spharmderivn(Ns, xi);

% Compute tangent basis.
d1 = bsxfun(@times, d1Y * cs, x) + bsxfun(@times, rho, d1s);
d2 = bsxfun(@times, d2Y * cs, x) + bsxfun(@times, rho, d2s);

end