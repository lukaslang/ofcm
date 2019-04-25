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
function tests = vbasisfunTest
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
h = 0;

% Create point x.
x = [0, 0, 1];
m = size(x, 1);

% Create evaluation points.
[~, y] = sphTriang(4);
n = size(y, 1);

% Create basis functions.
u = vbasisfun(k, h, x, y);
verifyEqual(testCase, size(u), [2*m, n, 3]);

end

function visualizeTest(testCase)

% Set parameters.
k = 3;
h = 0.75;

% Create point x.
x = [0, 0, 1];
m = size(x, 1);

% Create triangulation for visualization purpose.
[F, V] = sphTriang(4);
n = size(V, 1);

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create basis functions.
b = basisfun(k, h, x, V);
verifyEqual(testCase, size(b), [m, n]);

% Create basis functions.
y = vbasisfun(k, h, x, V);
verifyEqual(testCase, size(y), [2*m, n, 3]);

% Permute for plotting.
u = permute(y(1:m, :, :), [2, 3, 1]);
v = permute(y(m+1:end, :, :), [2, 3, 1]);

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), b, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
quiver3(V(:, 1), V(:, 2), V(:, 3), u(:, 1), u(:, 2), u(:, 3), 'g');
quiver3(V(:, 1), V(:, 2), V(:, 3), v(:, 1), v(:, 2), v(:, 3), 'w');

% Compute components of basis functions.
[dy1, dy2] = vbasiscompmem(k, h, x, V, 1024^3);

% Find coordinates.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
xi = [pi/2 - el, az];

% Compute coordinate basis at evaluation points.
[d1, d2] = sphtanbasis(xi, eye(3));
d1 = permute(d1, [3, 1, 2]);
d2 = permute(d2, [3, 1, 2]);

% Compute components of vectorial basis functions in canonical basis.
C1 = (dy1 .* d1(:, :, 1) + dy2 .* d2(:, :, 1))';
C2 = (dy1 .* d1(:, :, 2) + dy2 .* d2(:, :, 2))';
C3 = (dy1 .* d1(:, :, 3) + dy2 .* d2(:, :, 3))';

c = [1, 0]';
u = cat(3, C1 * c, C2 * c, C3 * c);
c = [0, 1]';
v = cat(3, C1 * c, C2 * c, C3 * c);

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), b, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
quiver3(V(:, 1), V(:, 2), V(:, 3), u(:, 1), u(:, 2), u(:, 3), 'g');
quiver3(V(:, 1), V(:, 2), V(:, 3), v(:, 1), v(:, 2), v(:, 3), 'w');

end