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
function [d11, d12, d21, d22] = sphtanbasisderiv(xi)
%SPHTANBASISDERIV Returns derivatives of tangent basis at evaluation points.
%   
%   [d11, d12, d21, d22] = SPHTANBASISDERIV(xi) takes spherical coordinates
%   xi and returns partial derivatives of tangent basis.
%
%   xi is a n-by-2 matrix with elements in [0, pi] and [0, 2*pi).
%   d11, d12, d21, d22 are n-by-3 matrices of row vectors on the 2-sphere.

n = size(xi, 1);
el = xi(:, 1);
az = xi(:, 2);  

d11 = - [sin(el).*cos(az), sin(el).*sin(az), cos(el)];
d12 = [-cos(el).*sin(az), cos(el).*cos(az), zeros(n, 1)];
d21 = [-cos(el).*sin(az), cos(el).*cos(az), zeros(n, 1)];
d22 = - [sin(el).*cos(az), sin(el).*sin(az), zeros(n, 1)];

end