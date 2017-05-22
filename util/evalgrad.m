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
function gradfx = evalgrad(f, scale, S, N, sc, bandwidth, layers)
%EVALGRAD Evaluates volumetric data on a given surface.
%   
%   gradfx = evalgrad(f, scale, S, sc, bandwidth, layers) takes a cell 
%   array f, scaling factors sc, evaluation points S, an offset sc, and 
%   bandwidth and layers specifying a narrow band, and returns gradfx at 
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
assert(all(size(scale) == [1, 3]));
assert(length(f) == length(S));
assert(length(f) == length(N));
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

% Compute gradient in embedding space.
[gy, gx, gz] = cellfun(@(x) gradient(double(permute(x, [2, 1, 3])), scale(1), scale(2), scale(3)), f, 'UniformOutput', false);

% Obtain data.
gx = cellfun(@(s, x) dataFromCube(s(:, 1), s(:, 2), s(:, 3), X, Y, Z, x), P, gx, 'UniformOutput', false);
gx = cellfun(@(x, s) mean(reshape(x, size(s, 1), length(rs)), 2), gx, S, 'UniformOutput', false);
gy = cellfun(@(s, x) dataFromCube(s(:, 1), s(:, 2), s(:, 3), X, Y, Z, x), P, gy, 'UniformOutput', false);
gy = cellfun(@(x, s) mean(reshape(x, size(s, 1), length(rs)), 2), gy, S, 'UniformOutput', false);
gz = cellfun(@(s, x) dataFromCube(s(:, 1), s(:, 2), s(:, 3), X, Y, Z, x), P, gz, 'UniformOutput', false);
gz = cellfun(@(x, s) mean(reshape(x, size(s, 1), length(rs)), 2), gz, S, 'UniformOutput', false);

% Compute dot product.
ip = cellfun(@(n, x, y, z) dot(n, [x, y, z], 2), N, gx, gy, gz, 'UniformOutput', false);

% Project to surface.
gradfx = cellfun(@(n, ipn, x, y, z) [x, y, z] - bsxfun(@times, ipn, n), N, ip, gx, gy, gz, 'UniformOutput', false);

end