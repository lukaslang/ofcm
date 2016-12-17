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
function [d1, d2] = sphtanbasis(xi, e)
%SPHTANBASIS Returns tangent basis at evaluation points.
%   
%   [d1, d2] = SPHTANBASIS(xi, e) takes spherical coordinates xi and an ONB
%   e of R^3 and returns tangent basis in the coordinate system specified
%   by e.
%
%   xi is a n-by-2 matrix with elements in [0, pi] and [0, 2*pi).
%   e is a 3-by-3 matrix where each row is a element of the basis.
%   d1, d2 are n-by-3 matrices of row vectors on the 2-sphere.

d1 = bsxfun(@times, cos(xi(:, 1)).*cos(xi(:, 2)), e(1, :));
d1 = d1 + bsxfun(@times, cos(xi(:, 1)).*sin(xi(:, 2)), e(2, :));
d1 = d1 + bsxfun(@times, -sin(xi(:, 1)), e(3, :));

d2 = bsxfun(@times, -sin(xi(:, 1)).*sin(xi(:, 2)), e(1, :));
d2 = d2 + bsxfun(@times, sin(xi(:, 1)).*cos(xi(:, 2)), e(2, :));

end