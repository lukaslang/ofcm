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

% This script executes all functionality of the software package.
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
path = fullfile('results', '2017-05-10-21-00-48');
file = fullfile(path, strcat(name, '.mat'));

% Load experiment.
load(file);

% Create triangulation for visualisation purpose.
[F, V] = halfsphTriang(7);
[S, ~] = cellfun(@(c) surfsynth(Ns, V, c), cs, 'UniformOutput', false);

% Evaluate data at vertices.
fd = evaldata(f, scale, S, sc, bandwidth, layers);

% Create segmentation.
sfd = cellfun(@(x) double(im2bw(x, graythresh(x))), fd, 'UniformOutput', false);
sfd = fd;
sfd = cellfun(@(x) ones(size(x, 1), 1), fd, 'UniformOutput', false);

% Find midpoints of faces on sphere.
TR = TriRep(F, V);
IC = normalise(TR.incenters);
[ICS, ~] = cellfun(@(c) surfsynth(Ns, IC, c), cs, 'UniformOutput', false);

% Evaluate basis functions at vertices.
[bfc1, bfc2] = vbasiscompmem(k, h, X, IC, mem);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(IC(:, 1), IC(:, 2), IC(:, 3));
el = pi/2 - el;

% Compute tangent basis.
[d1, d2] = cellfun(@(c) surftanbasis(Ns, c, [el, az]), cs(1:end-1), 'UniformOutput', false);

% Run through all pairs of frames.
for t=1:length(frames)-1
    
    % Compute pushforward of basis functions.
    v = bsxfun(@times, full((bfc1')*ofc{t}), d1{t}) + bsxfun(@times, full((bfc2')*ofc{t}), d2{t});

    figure(1);
    cla;
    hold on;
    colorbar;
    colormap(cmap);
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    daspect([1, 1, 1]);
    view(2);
    quiver3(ICS{t}(:, 1), ICS{t}(:, 2), ICS{t}(:, 3), v(:, 1), v(:, 2), v(:, 3), 0, 'r');
    title('Velocity field and surface data', 'FontName', 'Helvetica', 'FontSize', 14);

    figure(2);
    cla;
    hold on;
    colorbar;
    colormap(cmap);
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), sfd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    daspect([1, 1, 1]);
    view(2);
    title('Segmentation', 'FontName', 'Helvetica', 'FontSize', 14);

    figure(3);
    cla;
    hold on;
    colorbar;
    colormap(cmap);
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    daspect([1, 1, 1]);
    view(2);
    title('First frame.', 'FontName', 'Helvetica', 'FontSize', 14);

    figure(4);
    cla;
    hold on;
    colorbar;
    colormap(cmap);
    trisurf(F, S{t+1}(:, 1), S{t+1}(:, 2), S{t+1}(:, 3), fd{t+1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    daspect([1, 1, 1]);
    view(2);
    title('Second frame.', 'FontName', 'Helvetica', 'FontSize', 14);

    % Project and scale flow.
    up = projecttoplane(v);

    % Compute colour space scaling.
    nmax = max(sqrt(sum(up.^2, 2)));

    % Compute colour of projection.
    col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

    % Visualize flow.
    figure(5);
    cla;
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', col, 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(2);
    title('Colour-coded velocity.', 'FontName', 'Helvetica', 'FontSize', 14);

    % Plot colourwheel.
    figure(6);
    cla;
    cw = colourWheel;
    surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(2);
    title('Colour disk.', 'FontName', 'Helvetica', 'FontSize', 14);

    % Save figures.
    mkdir(fullfile(outputPath, name));
    export_fig(1, fullfile(outputPath, name, sprintf('%s-vel-frame-%.3i.png', name, t)), '-png', '-q300', '-a1', '-transparent');
    export_fig(2, fullfile(outputPath, name, sprintf('%s-seg-frame-%.3i.png', name, t)), '-png', '-q300', '-a1', '-transparent');
    export_fig(3, fullfile(outputPath, name, sprintf('%s-img1-frame-%.3i.png', name, t)), '-png', '-q300', '-a1', '-transparent');
    export_fig(4, fullfile(outputPath, name, sprintf('%s-img2-frame-%.3i.png', name, t)), '-png', '-q300', '-a1', '-transparent');
    export_fig(5, fullfile(outputPath, name, sprintf('%s-col-frame-%.3i.png', name, t)), '-png', '-q300', '-a1', '-transparent');
end