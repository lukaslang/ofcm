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

% This scrip renders a basis function.
clear;
close all;
clc;

% Define render quality (set to '-r600' for print quality).
quality = '-r300';

%% Zonal function.

% Set parameters of basis function.
k = 3;
h = 0.6;

% Create center point x of basis function.
x = [0, 0, 1];
m = size(x, 1);

% Create evaluation points.
[F, V] = sphTriang(5);
n = size(V, 1);

% Create basis functions.
b = basisfun(k, h, x, V);

figure(1);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), b, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
set(gca, 'ZLim', [-1, 1]);
set(gca, 'XLim', [-1, 1]);
set(gca, 'YLim', [-1, 1]);
set(gca, 'XTick', -1:0.5:1);
set(gca, 'YTick', -1:0.5:1);
set(gca, 'ZTick', -1:0.5:1);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);
colorbar;
title('Zonal function at north pole.');
export_fig(fullfile('results', 'basisfun.png'), '-png', quality, '-transparent', '-a1');

%% Vectorial basis.

% Create basis functions.
[Fv, Vv] = sphTriang(3);
y = vbasisfun(k, h, x, Vv);

% Permute for plotting.
u = permute(y(1:m, :, :), [2, 3, 1]);
v = permute(y(m+1:end, :, :), [2, 3, 1]);

figure(2);
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), b, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
set(gca, 'ZLim', [-1, 1]);
set(gca, 'XLim', [-1, 1]);
set(gca, 'YLim', [-1, 1]);
set(gca, 'XTick', -1:0.5:1);
set(gca, 'YTick', -1:0.5:1);
set(gca, 'ZTick', -1:0.5:1);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);
view(3);
title('Vectorial basis functions.');
colorbar;

quiver3(Vv(:, 1), Vv(:, 2), Vv(:, 3), u(:, 1), u(:, 2), u(:, 3), 'r');
quiver3(Vv(:, 1), Vv(:, 2), Vv(:, 3), v(:, 1), v(:, 2), v(:, 3), 'w');
export_fig(fullfile('results', 'vbasisfun.png'), '-png', quality, '-transparent', '-a1');