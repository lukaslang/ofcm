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
function [xi, w] = gausslegendre(deg, gamma)
%GAUSSLEGENDRE Returns evaluation points and weights for the Gauss-Legendre
%quadrature for a spherical cap of the 2-sphere.
%   
%   [xi, w] = GAUSSLEGENDRE(deg, gamma) takes an integer deg and a radius 
%   gamma in (0, pi), and returns evaluation points and weights for 
%   numerical integration of a spherical cap of radius gamma by means of
%   the Gauss-Legendre rule up to degree N.
%
%   xi and w are vectors of length N = (floor(deg/2) + 1)*(deg + 1).
assert(isscalar(deg));
assert(deg >= 0);
assert(isscalar(gamma));
assert(gamma < pi);
assert(gamma > 0);

% Compute m according to Gauss-Legendre.
m = floor(deg/2) + 1;

% Compute Gauss-Legendre quadrature rule.
[xi, w] = lgwt(m, -1, 1);

% Scale for spherical cap.
xi = (1 - cos(gamma)) * xi / 2 + (1 + cos(gamma)) / 2;
w = 2*pi*w*(1 - cos(gamma)) / (2*deg+2);

% Create tensor product quadrature rule.
xi1 = kron(acos(xi), ones(deg + 1, 1));

xi2 = 2*pi*(0:deg)/(deg + 1);
xi2 = kron(ones(m, 1), xi2');

% Create tensorial evaluation points.
xi = [xi1, mod(xi2, 2*pi)];

% Create tensorial weights.
w = kron(w, ones(deg + 1, 1));

end