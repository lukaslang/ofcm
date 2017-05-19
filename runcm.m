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
%frames = 139:141;
frames = 112:115;
%frames = 122:124;

% Load colormap.
load(fullfile('data', 'cmapblue.mat'));

% Specify max. memory for matrix multiplication.
mem = 3*1024^3;

% Define Gaussian filter.
sigma = 1.5;
hsize = 30;

% Define thresholding filter.
threshold = 80;
area = 100;

% Define surface fitting parameters.
Ns = 0:10;
beta0 = 1e-4;
beta1 = 10;
s = 3+eps;

% Set temporal derivative.
dt = 1;

% Parameters for radial projection of the data.
bandwidth = [0.9, 1.1];
layers = 40;

% Set subdivision parameter (number of basis functions is approx. 10*4^n).
ref = 5;

% Set parameters for basis function.
k = 3;
h = 0.99;

% Define degree of integration.
deg = 600;

% Set regularisation parameter.
alpha = 0.05;
beta = 0.001;
gamma = 0.005;

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
%s = cellfun(@(x) double(im2bw(x, graythresh(x))), fx, 'UniformOutput', false);
s = fx;
%s = cellfun(@(x) ones(size(x, 1), 1), fx, 'UniformOutput', false);

% Select frame.
t = 1;

% Compute optimality conditions.
[~, A, D, E, G, b] = optcondcm(Ns, cs{t}, cs{t+1}, X, k, h, xi, w, gradfx{t}, dtfx{t}, fx{t}, s{t}, mem);

% Solve linear system.
ofc = (A + alpha * D + beta * E + gamma * G) \ b;
fprintf('Relative residual %e.\n', norm((A + alpha * D + beta * E + gamma * G)*ofc - b)/norm(b));

% Create triangulation for visualization purpose.
[F, V] = halfsphTriang(6);
[S, ~] = cellfun(@(c) surfsynth(Ns, V, c), cs, 'UniformOutput', false);

% Evaluate data at vertices.
fd = evaldata(f, scale, S, sc, bandwidth, layers);

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

% Compute pushforward of basis functions.
v = bsxfun(@times, full((bfc1')*ofc), d1{t}) + bsxfun(@times, full((bfc2')*ofc), d2{t});

% Create segmentation.
%sfd = cellfun(@(x) double(im2bw(x, graythresh(x))), fd, 'UniformOutput', false);
sfd = fd;

figure;
hold on;
colormap(cmap);
trisurf(F, S{t+1}(:, 1), S{t+1}(:, 2), S{t+1}(:, 3), fd{t+1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(2);
quiver3(ICS{t+1}(:, 1), ICS{t+1}(:, 2), ICS{t+1}(:, 3), v(:, 1), v(:, 2), v(:, 3), 0, 'r');

figure;
hold on;
colormap(cmap);
trisurf(F, S{t+1}(:, 1), S{t+1}(:, 2), S{t+1}(:, 3), sfd{t+1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(2);

figure;
hold on;
colormap(cmap);
trisurf(F, S{t+1}(:, 1), S{t+1}(:, 2), S{t+1}(:, 3), fd{t+1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(2);
figure;
hold on;
colormap(cmap);
trisurf(F, S{t+1}(:, 1), S{t+1}(:, 2), S{t+1}(:, 3), fd{t+1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(2);

% Project and scale flow.
up = projecttoplane(v);

% Compute colour space scaling.
nmax = max(sqrt(sum(up.^2, 2)));

% Compute colour of projection.
col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

% Visualize flow.
figure;
hold on;
trisurf(F, S{t+1}(:, 1), S{t+1}(:, 2), S{t+1}(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', col, 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

% Plot colourwheel.
figure;
cw = colourWheel;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(2);