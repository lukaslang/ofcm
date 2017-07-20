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

% This scripts loads the data, fits a surface, extracts surface data, and
% plots the result.
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
path = fullfile(datapath, 'LSM 16.03.2012');
file = fullfile(path, strcat(name, '.lsm'));

% Load colormap.
load(fullfile('data', 'cmapblue.mat'));

% Select frames.
%frames = 112:115;
frames = 130:151;

% Define Gaussian filter.
sigma = 1.5;
hsize = 30;

% Define thresholding filter.
threshold = 80;
area = 100;

% Define surface fitting parameters.
Ns = 0:10;
beta0 = 1e-4;
beta1 = 100;
s = 3+eps;

% Parameters for radial projection of the data.
bandwidth = [0.9, 1.1];
layers = 40;

% Read dataset.
[f, scale] = loaddata(file, 1, frames);

% Scale.
scale = scale * 1e6;

% Reverse z-coordinate.
f = cellfun(@(X) flip(X, 3), f, 'UniformOutput', false);

% Find approximate cell centers.
CC = cellcenters(f, sigma, hsize, threshold, area);

% Scale local maxima.
C = cellfun(@(X) bsxfun(@times, X, scale), CC, 'UniformOutput', false);

% Align cell centers around origin.
[C, sc, sr] = aligncenters(C, 300);

% Fit surface and obtain Fourier coefficients.
[cs, ~] = fitsurface(Ns, C, beta0, beta1, s);

% Create triangulation for visualisation purpose.
[F, V] = halfsphTriang(7);
[S, rho] = cellfun(@(c) surfsynth(Ns, V, c), cs, 'UniformOutput', false);

% Normalise data.
f = cellfun(@(x) double(x) ./ 255, f, 'UniformOutput', false);

% Evaluate data.
fd = evaldata(f, scale, S, sc, bandwidth, layers);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;
xi = [el, az];

% Compute surface normals.
N = cellfun(@(x) surfnormals(Ns, x, xi), cs, 'UniformOutput', false);

% Compute surface gradient.
gradfd = evalgrad(f, scale, S, N, sc, bandwidth, layers);

% Select frame.
t = 1;

% Cell centres with fitted surface.
figure(1);
colormap(cmap);
hold on;
scatter3(C{t}(:, 1), C{t}(:, 2), C{t}(:, 3), '*r');
trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
axis square;
daspect([1, 1, 1]);
set(gca, 'ZLim', [0, 450]);
set(gca, 'XLim', [-450, 450]);
set(gca, 'YLim', [-450, 450]);
set(gca, 'XTick', -450:150:450);
set(gca, 'YTick', -450:150:450);
set(gca, 'ZTick', -450:150:450);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);
title('Cell centres with fitted surface and data.');
view(3);

% Raw data with surface.
figure(2);
colormap(cmap);
hold on;
vol3d('cdata', f{t}, 'XData', scale(1) * [0, 512] - sc(1), 'YData', scale(2) * [0, 512] - sc(2), 'ZData', scale(3) * [0, 44] - sc(3));
trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
axis square;
daspect([1, 1, 1]);
set(gca, 'ZLim', [0, 450]);
set(gca, 'XLim', [-450, 450]);
set(gca, 'YLim', [-450, 450]);
set(gca, 'XTick', -450:150:450);
set(gca, 'YTick', -450:150:450);
set(gca, 'ZTick', -450:150:450);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);
title('Raw data with surface and surface data.');
view(3);

% Gradient of surface data.
figure(3);
colormap(cmap);
hold on;
trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
quiver3(S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), gradfd{t}(:, 1), gradfd{t}(:, 2), gradfd{t}(:, 3), 1, 'r');
axis square;
daspect([1, 1, 1]);
set(gca, 'ZLim', [0, 450]);
set(gca, 'XLim', [-450, 450]);
set(gca, 'YLim', [-450, 450]);
set(gca, 'XTick', -450:150:450);
set(gca, 'YTick', -450:150:450);
set(gca, 'ZTick', -450:150:450);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);
title('Surface gradient of surface data.');
view(2);

% Animated surfac data.
figure(4);
colormap(cmap);
hold on;
axis square;
daspect([1, 1, 1]);
set(gca, 'ZLim', [0, 450]);
set(gca, 'XLim', [-450, 450]);
set(gca, 'YLim', [-450, 450]);
set(gca, 'XTick', -450:150:450);
set(gca, 'YTick', -450:150:450);
set(gca, 'ZTick', -450:150:450);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);
title('Raw data with surface and surface data.');
view(3);
for t=1:length(frames)
    cla;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    drawnow();
end