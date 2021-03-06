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
function tests = vbasiscompTest
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
k = 4;
h = 0.9;
ref = 3;
deg = 100;

% Create integration points and quadrature rule for spherical cap.
[xi, ~] = gausslegendre(deg, pi/2);
n = size(xi, 1);
Y = sphcoord(xi, eye(3));

% Triangulate of the upper unit hemi-sphere for placement of basis functions.
[~, X] = halfsphTriang(ref);
m = size(X, 1);

% Create basis functions.
tic;
[dy1, dy2] = vbasiscomp(k, h, X, Y);
toc;
verifyEqual(testCase, size(dy1), [2*m, n]);
verifyEqual(testCase, size(dy2), [2*m, n]);

% Create basis functions.
tic;
[dy1m, dy2m] = vbasiscompmem(k, h, X, Y, 128*8*m);
toc;
verifyEqual(testCase, size(dy1), [2*m, n]);
verifyEqual(testCase, size(dy2), [2*m, n]);

verifyEqual(testCase, dy1m, dy1, 'absTol', 1e-13);
verifyEqual(testCase, dy2m, dy2, 'absTol', 1e-13);

end

function visualizeTest(testCase)

% Set parameters.
k = 5;
h = 0.99;

% Create point x.
X = [0, 0, 1];
m = size(X, 1);

% Create evaluation points.
[F, Y] = sphTriang(5);
n = size(Y, 1);

% Find points at north and south pole.
idx = Y(:, 3) > 1-eps | Y(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
Y(idx, :) = Y(idx, :) * R;

[dy1, dy2] = vbasiscomp(k, h, X, Y);
dy1 = full(dy1);
dy2 = full(dy2);

% Create basis functions.
y = vbasisfun(k, h, X, Y);
verifyEqual(testCase, size(y), [2*m, n, 3]);

% Permute for plotting.
u = permute(y(1:m, :, :), [2, 3, 1]);
v = permute(y(m+1:end, :, :), [2, 3, 1]);

figure;
subplot(2, 2, 1);
hold on;
trisurf(F, Y(:, 1), Y(:, 2), Y(:, 3), dy1(1, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
quiver3(Y(:, 1), Y(:, 2), Y(:, 3), u(:, 1), u(:, 2), u(:, 3), 'g');

subplot(2, 2, 2);
hold on;
trisurf(F, Y(:, 1), Y(:, 2), Y(:, 3), dy2(1, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
quiver3(Y(:, 1), Y(:, 2), Y(:, 3), u(:, 1), u(:, 2), u(:, 3), 'g');

subplot(2, 2, 3);
hold on;
trisurf(F, Y(:, 1), Y(:, 2), Y(:, 3), dy1(2, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
quiver3(Y(:, 1), Y(:, 2), Y(:, 3), v(:, 1), v(:, 2), v(:, 3), 'w');

subplot(2, 2, 4);
hold on;
trisurf(F, Y(:, 1), Y(:, 2), Y(:, 3), dy2(2, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
quiver3(Y(:, 1), Y(:, 2), Y(:, 3), v(:, 1), v(:, 2), v(:, 3), 'w');

end