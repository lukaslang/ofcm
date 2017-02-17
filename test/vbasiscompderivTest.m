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
function tests = vbasiscompderivTest
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
tic;
[dy1d1, dy1d2, dy2d1, dy2d2] = vbasiscompderiv(k, h, x, y);
toc;
verifyEqual(testCase, size(dy1d1), [2*m, n]);
verifyEqual(testCase, size(dy1d2), [2*m, n]);
verifyEqual(testCase, size(dy2d1), [2*m, n]);
verifyEqual(testCase, size(dy2d2), [2*m, n]);

tic;
[dy1d1m, dy1d2m, dy2d1m, dy2d2m] = vbasiscompderivmem(k, h, x, y, 256*3*8*m);
toc;
verifyEqual(testCase, size(dy1d1m), [2*m, n]);
verifyEqual(testCase, size(dy1d2m), [2*m, n]);
verifyEqual(testCase, size(dy2d1m), [2*m, n]);
verifyEqual(testCase, size(dy2d2m), [2*m, n]);

verifyEqual(testCase, dy1d1, dy1d1m, 'absTol', 1e-16);
verifyEqual(testCase, dy1d2, dy1d2m, 'absTol', 1e-16);
verifyEqual(testCase, dy2d1, dy2d1m, 'absTol', 1e-16);
verifyEqual(testCase, dy2d2, dy2d2m, 'absTol', 1e-16);

end

function visualizeTest(testCase)

% Set parameters.
k = 4;
h = 0.75;

% Create point x.
x = [0, -1, 0];
m = size(x, 1);

% Create evaluation points.
[F, V] = sphTriang(5);
n = size(V, 1);

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
view(3);

% Compute partial derivatives.
[dy1d1, dy1d2, dy2d1, dy2d2] = vbasiscompderiv(k, h, x, V);
dy1d1 = full(dy1d1);
dy1d2 = full(dy1d2);
dy2d1 = full(dy2d1);
dy2d2 = full(dy2d2);

figure;
hold on;
subplot(1, 8, 1);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dy1d1(1, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
subplot(1, 8, 2);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dy1d2(1, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
subplot(1, 8, 3);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dy2d1(1, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
subplot(1, 8, 4);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dy2d2(1, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
subplot(1, 8, 5);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dy1d1(2, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
subplot(1, 8, 6);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dy1d2(2, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
subplot(1, 8, 7);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dy2d1(2, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
subplot(1, 8, 8);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), dy2d2(2, :), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);

end