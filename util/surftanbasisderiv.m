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
function [d11, d12, d21, d22] = surftanbasisderiv(Ns, cs, xi)
%SURFTANBASISDERIV Computes derivative of tangent basis on surface.
%   
%   [d11, d12, d21, d22] = SURFTANBASISDERIV(ns, cs, xi) takes degrees Ns 
%   and corresponding Fourier coefficients cs and evaluation points xi and 
%   returns the tangent basis.
%
%   Ns is a list os non-negative consecutive integers.
%   cs is a vector of Fourier coefficients.
%   xi is an n-by-2 matrix of coordinates [el, az], where el in [0, pi] and
%   az in [0, 2pi).
%
%   d11, d12, d21, d22 are n-by-3 matrices.

% Compute tangent basis of sphere.
[d1s, d2s] = sphtanbasis(xi, eye(3));

% Compute derivatives of tangent basis of sphere.
[d11s, ~, d21s, d22s] = sphtanbasisderiv(xi);

% Compute points on sphere.
x = sphcoord(xi, eye(3));

% Compute synthesis of evaluation points.
[~, rho] = surfsynth(Ns, x, cs);

% Compute partial derivatives of spherical harmonics.
[d1Y, d2Y] = spharmderivn(Ns, xi);

% Compute second partial derivatives of spherical harmonics.
[d11Y, ~, d21Y, d22Y] = spharmderivn2(Ns, xi);

% Compute partial derivatives of tangent basis.
d11 = bsxfun(@times, d11Y * cs, x) + 2 * bsxfun(@times, d1Y * cs, d1s) + bsxfun(@times, rho, d11s);
d21 = bsxfun(@times, d21Y * cs, x) + bsxfun(@times, d1Y * cs, d2s) + bsxfun(@times, d2Y * cs, d1s) + bsxfun(@times, rho, d21s);
d12 = d21;
d22 = bsxfun(@times, d22Y * cs, x) + 2 * bsxfun(@times, d2Y * cs, d2s) + bsxfun(@times, rho, d22s);

end