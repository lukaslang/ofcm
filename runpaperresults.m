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

% This script loads the data, computes optimality conditions, and outputs
% several files:
%
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]-data.mat
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]-linsys-000.mat
%
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
path = fullfile(datapath, 'LSM 16.03.2012');
file = fullfile(path, strcat(name, '.lsm'));

% Select frames.
frames = 50:151;

% Create start date and time.
startdate = datestr(now, 'yyyy-mm-dd-HH-MM-SS');

% Define and create output folder.
outputPath = fullfile('results', name);
mkdir(outputPath);

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
beta1 = 100;
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

% Normalise data.
f = cellfun(@(x) double(x) ./ 255, f, 'UniformOutput', false);

% Create segmentation function.
segfh = @(x) x;

% Plot 2D colourwheel.
figure(1);
cw = colourwheelbg;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(2);
%export_fig(fullfile(outputPath, 'colourwheel.png'), '-png', quality, '-a1', '-transparent');
%close all;

%% Optical flow vs. mass conservation.

% Select pair of frames.
t = 62;

% Compute synthesis of evaluation points.
Sy1 = surfsynth(Ns, Y, cs{t});
Sy2 = surfsynth(Ns, Y, cs{t+1});

% Evaluate data.
fx1 = cell2mat(evaldata(f(t), scale, {Sy1}, sc, bandwidth, layers));
fx2 = cell2mat(evaldata(f(t+1), scale, {Sy2}, sc, bandwidth, layers));

% Compute temporal derivative.
dtfx = (fx2 - fx1) / dt;
clear fx2;

% Compute surface normals.
N = surfnormals(Ns, cs{t}, xi);

% Compute surface gradient.
gradfx = cell2mat(evalgrad(f(t), scale, {Sy1}, {N}, sc, bandwidth, layers));
clear N;
clear Sy1;
clear Sy2;

fprintf('Computing optimaly conditions.\n');
timerVal = tic;
[~, Aof, Acm, D, E, G, bof, bcm] = optcondofcm(Ns, cs{t}, cs{t+1}, X, k, h, xi, w, gradfx, dtfx, fx1, segfh(fx1), mem);
elapsed = toc(timerVal);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Set regularisation parameters.
alpha = 0.01;
beta = 0.001;

% Solve linear system for optical flow.
timerVal = tic;
c = (Aof + alpha * D + beta * E) \ bof;
relres = norm((Aof + alpha * D + beta * E) * c - bof) / norm(bof);
elapsed = toc(timerVal);
fprintf('Relative residual %e.\n', relres);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Compute surface normals.
N = surfnormals(Ns, cs{t}, xi);

% Compute surface velocity.
dtrho = surfsynth(Ns, IC, (cs{t+1} - cs{t}) / dt);
Vs = bsxfun(@times, dtrho, IC);

% Compute flow.
w = bsxfun(@times, full((bfc1')*c{p}), d1) + bsxfun(@times, full((bfc2')*c{p}), d2);

% Project and scale flow.
wp = projecttoplane(w);

% Recover total velocity.
U = Vs + w;

% Project and scale flow.
Up = projecttoplane(U);

% Compute colour of projection.
col = double(squeeze(computeColour(wp(:, 1), wp(:, 2)))) ./ 255;

% Plot colour-coded flow.
figure(1);
cla;
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
view(3);
adjust3dplot;
%export_fig(fullfile(outputPath, 'of-flow3', folderstr, sprintf('flow3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
delete(C);

% Top view.
view(2);
C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
%export_fig(fullfile(outputPath, 'of-flow2', folderstr, sprintf('flow2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');

% Compute colour of projection.
col = double(squeeze(computeColour(Up(:, 1), Up(:, 2)))) ./ 255;

% Total velocity.
figure(1);
cla;
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
view(3);
adjust3dplot;
%export_fig(fullfile(outputPath, 'of-motion3', folderstr, sprintf('motion3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
delete(C);

% Top view.
view(2);
C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
%export_fig(fullfile(outputPath, 'of-motion2', folderstr, sprintf('motion2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');

% Visualise closeup.

% Set regularisation parameters.
alpha = 0.01;
beta = 0.001;
gamma = 0.005;

% Solve linear system for mass conservation.
timerVal = tic;
c = (Acm + alpha * D + beta * E + gamma * G) \ bcm;
relres = norm((Acm + alpha * D + beta * E + gamma * G) * c - bcm) / norm(bcm);
elapsed = toc(timerVal);
fprintf('Relative residual %e.\n', relres);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Visualise differences for different segmentations.

% Create segmentation function.
%segfh = @(x) double(im2bw(x, graythresh(x)));
segfh = @(x) x;
%segfh = @(x) ones(size(x, 1), 1);

