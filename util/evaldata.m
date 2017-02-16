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
function fx = evaldata(f, scale, S, sc, bandwidth, layers)
%EVALDATA Evaluates volumetric data on a given surface.
%   
%   fx = evaldata(f, scale, S, sc, bandwidth, layers) takes a cell 
%   array f, scaling factors sc, evaluation points S, an offset sc, and 
%   bandwidth and layers specifying a narrow band, and returns data fx at 
%   points S.
%
%   f is a cell array of length k with each containing a matrix [m, 3].
%   scale is a vector [sx, sy, sz].
%   S is a cell array of length k with each containing a matrix [n, 3].
%   sc is a vector [scx, scy, scz].
%   dt is a scalar > 0.
%   bandwidth is a vector [bandwidth(1), bandwidth(2)].
%   layers is a positive integer.
%
%   fx is a cell array of length k with each containing a vector.
assert(length(f) > 1);
assert(all(size(scale) == [1, 3]));
assert(length(f) == length(S));
assert(all(size(sc) == [1, 3]));

% Create grid.
[m, n, o] = size(f{1});
[X, Y, Z] = ndgrid(1:m, 1:n, 1:o);

% Scale grid.
X = scale(1) * X;
Y = scale(2) * Y;
Z = scale(3) * Z;

% Define narrow band around surface.
rs = linspace(bandwidth(1), bandwidth(2), layers);
P = cellfun(@(x) kron(rs', x), S, 'UniformOutput', false);

% Compute offset.
P = cellfun(@(s) bsxfun(@plus, s, sc), P, 'UniformOutput', false);

% Obtain data.
fx = cellfun(@(s, x) dataFromCube(s(:, 1), s(:, 2), s(:, 3), X, Y, Z, double(x)), P, f, 'UniformOutput', false);
fx = cellfun(@(x, s) max(reshape(x, size(s, 1), length(rs)), [], 2), fx, S, 'UniformOutput', false);

end