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
function test_suite = vbasiscompderivTest
    initTestSuite;
end

function resultTest

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
[dy1d1, dy1d2, dy2d1, dy2d2] = vbasiscompderiv(k, h, x, y);
assertEqual(size(dy1d1), [2*m, n]);
assertEqual(size(dy1d2), [2*m, n]);
assertEqual(size(dy2d1), [2*m, n]);
assertEqual(size(dy2d2), [2*m, n]);

end

function visualizeTest

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
assertEqual(size(b), [m, n]);

% Create basis functions.
y = vbasisfun(k, h, x, V);
assertEqual(size(y), [2*m, n, 3]);

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