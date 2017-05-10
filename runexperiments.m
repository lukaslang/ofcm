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
path = fullfile(datapath, 'LSM 16.03.2012');
file = fullfile(path, strcat(name, '.lsm'));

% Select frames.
frames = 112:115;

% Create start date and time.
startdate = datestr(now, 'yyyy-mm-dd-HH-MM-SS');

% Define and create output folder.
outputPath = fullfile('results', startdate);
mkdir(outputPath);

% Load colormap.
load(fullfile('data', 'cmapblue.mat'));

% Specify max. memory for matrix multiplication.
mem = 3*1024^3;

% Define Gaussian filter.
sigma = 1.5;
hsize = 30;

% Define thresholding filter.
t = 80;
area = 100;

% Define surface fitting parameters.
Ns = 0:10;
beta0 = 1e-4;
beta1 = 1;
s = 3+eps;

% Set temporal derivative.
dt = 1;

% Parameters for radial projection of the data.
bandwidth = [0.8, 1.2];
layers = 80;

% Set subdivision parameter (number of basis functions is approx. 10*4^n).
ref = 5;

% Set parameters for basis function.
k = 3;
h = 0.95;

% Define degree of integration.
deg = 400;

% Set regularisation parameter.
alpha = 0.01;
beta = 0.001;
gamma = 0.1;

% Read dataset.
[f, scale] = loaddata(file, 1, frames);

% Scale.
scale = scale * 1e6;

% Reverse z-coordinate.
f = cellfun(@(X) flip(X, 3), f, 'UniformOutput', false);

% Find approximate cell centers.
CC = cellcenters(f, sigma, hsize, t, area);

% Scale local maxima.
C = cellfun(@(X) bsxfun(@times, X, scale), CC, 'UniformOutput', false);

% Align cell centers around origin.
[C, sc, sr] = aligncenters(C, 300);

% Fit surface and obtain Fourier coefficients.
[cs, ~] = fitsurface(Ns, C, beta0, beta1, s);

% Triangulate of the upper unit hemi-sphere for placement of basis functions.
[~, X] = halfsphTriang(ref);

% Create integration points and quadrature rule for spherical cap.
[xi, w] = gausslegendre(deg, pi/2);

% Compute evaluation points.
Y = sphcoord(xi, eye(3));

% Compute synthesis of evaluation points.
[Sy, ~] = cellfun(@(c) surfsynth(Ns, Y, c), cs, 'UniformOutput', false);

% Normalise data.
f = cellfun(@(x) double(x) ./ 255, f, 'UniformOutput', false);

% Evaluate data.
fx = evaldata(f, scale, Sy, sc, bandwidth, layers);

% Compute temporal derivative.
dtfx = cellfun(@(x, y) (y - x) / dt, fx(1:end - 1), fx(2:end), 'UniformOutput', false);

% Compute surface normals.
N = cellfun(@(x) surfnormals(Ns, x, xi), cs, 'UniformOutput', false);

% Compute surface gradient.
gradfx = evalgrad(f, scale, Sy, N, sc, bandwidth, layers);

% Create segmentation.
s = cellfun(@(x) double(imbinarize(x)), fx, 'UniformOutput', false);
s = fx;

% Run through all pair of frames.
for t=1:length(frames)-1
    fprintf('Computing velocity field %i/%i.\n', t, length(frames)-1);
    
    % Compute optimality conditions.
    [~, A, D, E, G, b] = optcondcm(Ns, cs{t}, cs{t+1}, X, k, h, xi, w, gradfx{t}, dtfx{t}, fx{t}, s{t}, mem);

    % Solve linear system.
    [ofc{t}, L{t}] = solvesystem(A + alpha * D + beta * E + gamma * G, b, 1e-6, 2000);
    fprintf('GMRES terminated at iteration %i with relative residual %e.\n', L.iter(2), L.relres);

end

% Create triangulation for visualization purpose.
[F, V] = halfsphTriang(7);
[S, ~] = cellfun(@(c) surfsynth(Ns, V, c), cs, 'UniformOutput', false);

% Evaluate data at vertices.
fd = evaldata(f, scale, S, sc, bandwidth, layers);

% Create segmentation.
sfd = cellfun(@(x) double(imbinarize(x)), fd, 'UniformOutput', false);
sfd = fd;

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
    hold on;
    colorbar;
    colormap(cmap);
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    daspect([1, 1, 1]);
    view(2);
    quiver3(ICS{t}(:, 1), ICS{t}(:, 2), ICS{t}(:, 3), v(:, 1), v(:, 2), v(:, 3), 0, 'r');
    title('Velocity field and surface data', 'FontName', 'Helvetica', 'FontSize', 14);

    figure(2);
    hold on;
    colorbar;
    colormap(cmap);
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), sfd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    daspect([1, 1, 1]);
    view(2);
    title('Segmentation', 'FontName', 'Helvetica', 'FontSize', 14);

    figure(3);
    hold on;
    colorbar;
    colormap(cmap);
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    daspect([1, 1, 1]);
    view(2);
    title('First frame.', 'FontName', 'Helvetica', 'FontSize', 14);

    figure(4);
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
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', col, 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    title('Colour-coded velocity.', 'FontName', 'Helvetica', 'FontSize', 14);

    % Plot colourwheel.
    figure(6);
    cw = colourWheel;
    surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(2);
    title('Colour disk.', 'FontName', 'Helvetica', 'FontSize', 14);

    export_fig(1, fullfile(outputPath, name, sprintf('%s-frame-%.3i-vel.png', name, t)), '-png', '-q300', '-a1', '-transparent');
    export_fig(2, fullfile(outputPath, name, sprintf('%s-frame-%.3i-seg.png', name, t)), '-png', '-q300', '-a1', '-transparent');
    export_fig(3, fullfile(outputPath, name, sprintf('%s-frame-%.3i-img1.png', name, t)), '-png', '-q300', '-a1', '-transparent');
    export_fig(4, fullfile(outputPath, name, sprintf('%s-frame-%.3i-img2.png', name, t)), '-png', '-q300', '-a1', '-transparent');
    export_fig(5, fullfile(outputPath, name, sprintf('%s-frame-%.3i-col.png', name, t)), '-png', '-q300', '-a1', '-transparent');

end