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

% This script renders experiments and exports figures to the folder:
%
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]/[yyyy-mm-dd-HH-MM-SS]/
%
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';

% Set datestring of generated data.
timestamp1 = '2017-05-22-15-43-34';

% Set datestring of experiment.
timestamp2 = '2017-05-22-23-45-03';

% Define render quality.
quality = '-r300';

% Load data.
path = fullfile('results', name);
file = fullfile(path, sprintf('%s-data.mat', timestamp1));
load(file);

% Load colormap.
load(fullfile('data', 'cmapblue.mat'));

% Specify max. memory for matrix multiplication.
mem = 3*1024^3;

% Create triangulation for visualisation purpose.
[F, V] = halfsphTriang(5);
[S, rho] = cellfun(@(c) surfsynth(Ns, V, c), cs, 'UniformOutput', false);

% Evaluate data at vertices.
fd = evaldata(f, scale, S, sc, bandwidth, layers);

% Find midpoints of faces on sphere.
TR = TriRep(F, V);
IC = normalise(TR.incenters);
[ICS, ICrho] = cellfun(@(c) surfsynth(Ns, IC, c), cs, 'UniformOutput', false);

% Evaluate basis functions at vertices.
[bfc1, bfc2] = vbasiscompmem(k, h, X, IC, mem);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(IC(:, 1), IC(:, 2), IC(:, 3));
el = pi/2 - el;
xi = [el, az];

% Compute tangent basis.
[d1, d2] = cellfun(@(c) surftanbasis(Ns, c, xi), cs(1:end-1), 'UniformOutput', false);

% Compute tangent basis of sphere.
[d1s, d2s] = sphtanbasis(xi, eye(3));

% Precompute.
for t=1:length(frames)-1
    % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-of-%.3i.mat', timestamp2, frames(t)));
    load(file);
    
    for p=1:length(c)
        % Compute flow.
        v{t, p} = bsxfun(@times, full((bfc1')*c{p}), d1{t}) + bsxfun(@times, full((bfc2')*c{p}), d2{t});
        % Compute flow on shpere.
        vs{t, p} = bsxfun(@times, full((bfc1')*c{p}), d1s) + bsxfun(@times, full((bfc2')*c{p}), d2s);
    end
end

% Find min and max values.
for p=1:length(c)
    vmin{p} = min(min([v{:, p}]));
    vmax{p} = max(max([v{:, p}]));
    vsmin{p} = min(min([vs{:, p}]));
    vsmax{p} = max(max([vs{:, p}]));
    
    % Compute mean.
    meanvs{p} = mean(cat(3, vs{:, p}), 3);
end

% Select experiment.
p = 3;

% Project and scale.
up = projecttoplane(meanvs{p});

% Compute colour space scaling.
nmax = max(sqrt(sum(up.^2, 2)));

% Compute colour of projection.
col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

% Plot mean flow.
figure(1);
cla;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
view(2);
daspect([1, 1, 1]);
set(gca, 'ZLim', [0, 1]);
set(gca, 'XLim', [-1, 1]);
set(gca, 'YLim', [-1, 1]);
set(gca, 'XTick', -1:0.5:1);
set(gca, 'YTick', -1:0.5:1);
set(gca, 'ZTick', -1:0.5:1);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);
%export_fig(fullfile('results', 'of-meanflow.png'), '-png', quality, '-a1', '-transparent');

% Set seed points for streamlines.
[X, Y] = meshgrid(-1:0.05:1, -1:0.05:1);
idx = find(X.^2 + Y.^2 <= 1);
xs = [X(idx), Y(idx)];

% Set parameters.
nmax = max(sqrt(sum(meanvs{p}.^2, 2)));
h = 0.01/nmax;
maxit = 50;
lw = 1;

% Plot streamlines for meanflow.
figure(2);
cla;
hold on;
view(2);
daspect([1, 1, 1]);
set(gca, 'ZLim', [0, 1]);
set(gca, 'XLim', [-1, 1]);
set(gca, 'YLim', [-1, 1]);
set(gca, 'XTick', -1:0.5:1);
set(gca, 'YTick', -1:0.5:1);
set(gca, 'ZTick', -1:0.5:1);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);
colormap('summer');
streamlines2(IC, meanvs{p}, xs, h, maxit, 'summer', lw);
%export_fig(fullfile('results', 'of-streamlines.png'), '-png', quality, '-a1', '-transparent');

% Select frame.
t = length(frames)-1;

% Project and scale.
up = projecttoplane(v{t, p});

% Compute colour space scaling.
nmax = max(sqrt(sum(up.^2, 2)));

% Compute colour of projection.
col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

% Compute angles.
[theta, ~] = cart2pol(up(:, 1), up(:, 2));

% Plot rose diagram.
figure(3);
subplot(1, 2, 1);
cla;
hold on;
trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
adjust3dplot;
view(2);
subplot(1, 2, 2);
polarhistogram(theta, 20, 'Normalization', 'probability');