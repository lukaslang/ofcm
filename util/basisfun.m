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
function b = basisfun(k, h, x, y)
%BASISFUN Evaluates basis functions.
%   
%   b = BASISFUN(k, h, x, y) takes a degree k, a parameter h, a matrix 
%   x of locations of the basis functions, and a list of evaluation 
%   points y, and returns basis functions evaluated at y.
%
%   k is a non-negative integer, h a scalar in (0, 1).
%   x is a m-by-3 matrix where each row is a point on the 2-sphere.
%   y is a n-by-3 matrix where each row is a point on the 2-sphere.

% Compute dot products between all points.
b = x*(y');

% Evaluate basis functions at points xi.
idx = b > h;
b(~idx) = 0;
b(idx) = ((b(idx) - h) .^ k) / (1 - h)^k;

end