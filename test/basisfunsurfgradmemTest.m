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
function tests = basisfunsurfgradmemTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Set parameters.
k = 3;
h = 0.5;

% Create point x.
[~, x] = sphTriang(2);
m = size(x, 1);

% Create evaluation points.
[~, y] = sphTriang(4);
n = size(y, 1);

% Specify max. memory for matrix multiplication.
mem = 1*1024^3;

% Create basis functions.
bmem = basisfunsurfgradmem(k, h, x, y, mem);
assertEqual(testCase, size(bmem), [m, n]);

end

function visualizeTest(testCase)

% Set parameters.
k = 3;
h = 0.8;

% Create point where basis function is located.
X = [0, 0, 1];

% Create triangulation for visualization purpose.
[F, V] = sphTriang(4);

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Specify max. memory for matrix multiplication.
mem = 1024^3;

% Compute basis functions.
b = basisfunmem(k, h, X, V, mem);

% Create surface gradient of basis functions.
[b1, b2] = basisfunsurfgradmem(k, h, X, V, mem);

% Select a few of the vectors randomly for plotting.
rs = rand(size(V, 1), 1);
rs = find(rs <= 1);

% Compute coordinate basis at evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;
[d1, d2] = sphtanbasis([el, az], eye(3));

% Compute and plot surface gradient.
v = bsxfun(@times, b1', d1) + bsxfun(@times, b2', d2);

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), full(b), 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;
quiver3(V(rs, 1), V(rs, 2), V(rs, 3), v(rs, 1), v(rs, 2), v(rs, 3), 1, 'r', 'LineWidth', 2);
daspect([1, 1, 1]);

end