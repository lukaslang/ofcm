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
function x = sphcoord(xi, e)
%SPHCOORD Returns points in a spherical coordiante system.
%   
%   x = SPHCOORD(xi, e) takes spherical coordinates xi and an ONB e of R^3
%   and returns points on the 2-sphere in the coordinate system specified
%   by e.
%
%   xi is a n-by-2 matrix with elements in [0, pi] and [0, 2*pi).
%   e is a 3-by-3 matrix where each row is a element of the basis.
%   x is a n-by-3 matrix of row vectors on the 2-sphere.

% Compute evaluation points in R^3.
x = bsxfun(@times, sin(xi(:, 1)).*cos(xi(:, 2)), e(1, :));
x = x + bsxfun(@times, sin(xi(:, 1)).*sin(xi(:, 2)), e(2, :));
x = x + bsxfun(@times, cos(xi(:, 1)), e(3, :));

end