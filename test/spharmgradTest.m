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
function tests = spharmgradTest
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
N = 3;

% Create evaluation points.
[~, y] = sphTriang(4);
n = size(y, 1);

% Create basis functions.
u = spharmgrad(N, y);
verifyEqual(testCase, size(u), [n, 2*N + 1, 3]);

end

function visualizeTest(testCase)

% Set parameters.
N = 3;

% Create evaluation points.
[F, V] = sphTriang(4);
n = size(V, 1);

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create basis functions.
Ynj = spharm(N, V);

% Create basis functions.
u = spharmgrad(N, V);
verifyEqual(testCase, size(u), [n, 2*N + 1, 3]);

% Pick one order for plotting.
k = 2;
u = squeeze(u(:, k, :));

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), Ynj(:, k), 'EdgeColor', 'none', 'FaceColor', 'interp');
quiver3(V(:, 1), V(:, 2), V(:, 3), u(:, 1), u(:, 2), u(:, 3), 1, 'g');
daspect([1, 1, 1]);
view(3);

end