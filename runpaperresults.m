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

% Define render quality (set to '-r600' for print quality).
quality = '-r100';

% Create start date and time.
startdate = datestr(now, 'yyyy-mm-dd-HH-MM-SS');

% Define and create output folder.
outputPath = fullfile('results', name, 'figures', startdate);
mkdir(outputPath);

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
beta1 = 100;
s = 3+eps;

% Set temporal interval.
dt = 1;

% Parameters for radial projection of the data.
bandwidth = [0.9, 1.1];
layers = 40;

% Set subdivision parameter (number of basis functions is approx. 10*4^n).
ref = 5;
% ref = 6;

% Set parameters for basis function.
k = 3;
h = 0.995;

% Define degree of integration.
deg = 400;

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
[xi, weights] = gausslegendre(deg, pi/2);

% Compute evaluation points.
Y = sphcoord(xi, eye(3));

% Normalise data.
f = cellfun(@(x) double(x) ./ 255, f, 'UniformOutput', false);

% Create segmentation function.
segfh = @(x) x;

% Create triangulation for visualization purpose.
[F, V] = halfsphTriang(7);

% Find midpoints of faces on sphere.
TR = TriRep(F, V);
IC = normalise(TR.incenters);

%% Plot dataset and fitted surface.

% Create spherical mesh for visualisation and remove poles.
nsphere = 31;
[Xm, Ym, Zm] = sphere(nsphere-1);
Xm = Xm(2:end-1, :);
Ym = Ym(2:end-1, :);
Zm = Zm(2:end-1, :);

% Define scaling factor of spherical mesh.
gridscale = 1.01;

% Compute synthesis at mesh vertices.
[Vm, ~] = cellfun(@(c) surfsynth(Ns, [Xm(:), Ym(:), Zm(:)], c), cs, 'UniformOutput', false);
Vm = cellfun(@(x) x .* gridscale, Vm, 'UniformOutput', false);
Xm = cellfun(@(x) reshape(x(:, 1), nsphere-2, nsphere), Vm, 'UniformOutput', false);
Ym = cellfun(@(x) reshape(x(:, 2), nsphere-2, nsphere), Vm, 'UniformOutput', false);
Zm = cellfun(@(x) reshape(x(:, 3), nsphere-2, nsphere), Vm, 'UniformOutput', false);

% Select frames to render.
selFrames = 1:50:length(frames);

for t=selFrames
    % Compute synthesis for vertices.
    [S, rho] = surfsynth(Ns, V, cs{t});

    % Evaluate data.
    fd = cell2mat(evaldata(f(t), scale, {S}, sc, bandwidth, layers));

    % Plot data.
    plotdata(F, V, S, Xm{t}, Ym{t}, Zm{t}, f{t}, fd, rho, cmap, scale, sc, bandwidth);
    
    % Save figures.
    export_fig(fullfile(outputPath, sprintf('rho2-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1', figure(1));
    export_fig(fullfile(outputPath, sprintf('raw3-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1', figure(2));
    export_fig(fullfile(outputPath, sprintf('raw3-%s-surf-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1', figure(3));
    export_fig(fullfile(outputPath, sprintf('raw3-%s-surf-grid-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1', figure(4));
    export_fig(fullfile(outputPath, sprintf('raw3-%s-cross-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1', figure(5));
    export_fig(fullfile(outputPath, sprintf('data3-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1', figure(6));
    export_fig(fullfile(outputPath, sprintf('data2-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1', figure(7));
    export_fig(fullfile(outputPath, sprintf('data3-%s-grid-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1', figure(8));
    export_fig(fullfile(outputPath, sprintf('data2-%s-grid-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1', figure(9));
    close all;
end
clear Vm;
clear Xm;
clear Ym;
clear Zm;

%% Compute optimality conditions.

% Select pair of frames.
t = 63;

% Select data.
f1 = f(t);
f2 = f(t+1);
clear f;

% Compute synthesis of evaluation points.
Sy1 = surfsynth(Ns, Y, cs{t});
Sy2 = surfsynth(Ns, Y, cs{t+1});

% Evaluate data.
fx1 = cell2mat(evaldata(f1, scale, {Sy1}, sc, bandwidth, layers));
fx2 = cell2mat(evaldata(f2, scale, {Sy2}, sc, bandwidth, layers));

% Compute temporal derivative.
dtfx = (fx2 - fx1) / dt;
clear fx2;

% Compute surface normals.
N = surfnormals(Ns, cs{t}, xi);

% Compute surface gradient.
gradfx = cell2mat(evalgrad(f1, scale, {Sy1}, {N}, sc, bandwidth, layers));
clear N;
clear Sy1;
clear Sy2;

fprintf('Computing optimaly conditions.\n');
timerVal = tic;
[~, Aof, Acm, D, E, G, bof, bcm] = optcondofcm(Ns, cs{t}, cs{t+1}, X, k, h, xi, weights, gradfx, dtfx, fx1, segfh(fx1), mem);
elapsed = toc(timerVal);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

%% Create triangulation etc. for visualisation.

% Compute synthesis for vertices.
[S, ~] = surfsynth(Ns, V, cs{t});

% Evaluate data at vertices.
fd1 = cell2mat(evaldata(f1, scale, {S}, sc, bandwidth, layers));
fd2 = cell2mat(evaldata(f2, scale, {S}, sc, bandwidth, layers));

% Compute synthesis for midpoints.
[ICS, ~] = surfsynth(Ns, IC, cs{t});

% Compute surface normals.
N = surfnormals(Ns, cs{t}, xi);

% Compute surface velocity.
dtrho = surfsynth(Ns, IC, (cs{t+1} - cs{t}) / dt);
Vs = bsxfun(@times, dtrho, IC);

% Evaluate basis functions at vertices.
[bfc1, bfc2] = vbasiscompmem(k, h, X, IC, mem);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(IC(:, 1), IC(:, 2), IC(:, 3));
el = pi/2 - el;
xi = [el, az];

% Compute tangent basis.
[d1, d2] = surftanbasis(Ns, cs{t}, xi);

%% Compute optical flow.

% Set regularisation parameters.
alpha = 0.1;
beta = 0.001;

% Create folder string with parameters.
folderstr = sprintf('alpha-%.4g-beta-%.4g', alpha, beta);

% Solve linear system for optical flow.
timerVal = tic;
c = (Aof + alpha * D + beta * E) \ bof;
relres = norm((Aof + alpha * D + beta * E) * c - bof) / norm(bof);
elapsed = toc(timerVal);
fprintf('Relative residual %e.\n', relres);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Compute flow.
w = bsxfun(@times, full((bfc1')*c), d1) + bsxfun(@times, full((bfc2')*c), d2);

% Plot flow.
plotflow(F, S, ICS, fd1, fd2, w, cmap);

% Save figures.
export_fig(fullfile(outputPath, sprintf('of-flow3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(1));
export_fig(fullfile(outputPath, sprintf('of-flow2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(2));
export_fig(fullfile(outputPath, sprintf('of-flow2-detail1-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(3));
export_fig(fullfile(outputPath, sprintf('of-flow2-detail1-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(4));
export_fig(fullfile(outputPath, sprintf('of-flow2-detail2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(5));
export_fig(fullfile(outputPath, sprintf('of-flow2-detail2-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(6));
close all;

%% Compute total motion.

% Recover total velocity.
U = Vs + w;

% Plot flow.
plotflow(F, S, ICS, fd1, fd2, U, cmap);

% Save figures.
export_fig(fullfile(outputPath, sprintf('of-motion3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(1));
export_fig(fullfile(outputPath, sprintf('of-motion2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(2));
close all;

%% Mass conservation.

% Set regularisation parameters.
alpha = 0.01;
beta = 0.001;
gamma = 0.005;

% Create folder string with parameters.
folderstr = sprintf('alpha-%.4g-beta-%.4g-gamma-%.4g', alpha, beta, gamma);

% Solve linear system for mass conservation.
timerVal = tic;
c = (Acm + alpha * D + beta * E + gamma * G) \ bcm;
relres = norm((Acm + alpha * D + beta * E + gamma * G) * c - bcm) / norm(bcm);
elapsed = toc(timerVal);
fprintf('Relative residual %e.\n', relres);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Compute flow.
u = bsxfun(@times, full((bfc1')*c), d1) + bsxfun(@times, full((bfc2')*c), d2);

% Plot flow.
plotflow(F, S, ICS, fd1, fd2, u, cmap);

% Save figures.
export_fig(fullfile(outputPath, sprintf('cm-flow3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(1));
export_fig(fullfile(outputPath, sprintf('cm-flow2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(2));
export_fig(fullfile(outputPath, sprintf('cm-flow2-detail1-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(3));
export_fig(fullfile(outputPath, sprintf('cm-flow2-detail1-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(4));
export_fig(fullfile(outputPath, sprintf('cm-flow2-detail2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(5));
export_fig(fullfile(outputPath, sprintf('cm-flow2-detail2-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(6));
close all;

%% Compute total motion.

% Compute surface velocity.
dtrho = surfsynth(Ns, IC, (cs{t+1} - cs{t}) / dt);
Vs = bsxfun(@times, dtrho, IC);

% Compute scalar normal part of surface velocity.
Vsn = dot(Vs, N, 2);

% Recover total velocity.
U = bsxfun(@times, Vsn, N) + u;

% Plot flow.
plotflow(F, S, ICS, fd1, fd2, U, cmap);

% Save figures.
export_fig(fullfile(outputPath, sprintf('cm-motion3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(1));
export_fig(fullfile(outputPath, sprintf('cm-motion2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(2));
close all;

%% Use automatic thresholding for segmentation.

% Use automatic thresholding.
segfh = @(x) double(imbinarize(x, graythresh(x)));

fprintf('Computing optimaly conditions.\n');
timerVal = tic;
[~, Aof, Acm, D, E, G, bof, bcm] = optcondofcm(Ns, cs{t}, cs{t+1}, X, k, h, xi, weights, gradfx, dtfx, fx1, segfh(fx1), mem);
elapsed = toc(timerVal);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Set regularisation parameters.
alpha = 0.1;
beta = 0.001;

% Create folder string with parameters.
folderstr = sprintf('alpha-%.4g-beta-%.4g', alpha, beta);

% Solve linear system for optical flow.
timerVal = tic;
c = (Aof + alpha * D + beta * E) \ bof;
relres = norm((Aof + alpha * D + beta * E) * c - bof) / norm(bof);
elapsed = toc(timerVal);
fprintf('Relative residual %e.\n', relres);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Compute flow.
w = bsxfun(@times, full((bfc1')*c), d1) + bsxfun(@times, full((bfc2')*c), d2);

% Plot flow.
plotflow(F, S, ICS, fd1, fd2, w, cmap);

% Save figures.
export_fig(fullfile(outputPath, sprintf('of-flow3-graythresh-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(1));
export_fig(fullfile(outputPath, sprintf('of-flow2-graythresh-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(2));
export_fig(fullfile(outputPath, sprintf('of-flow2-graythresh-detail1-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(3));
export_fig(fullfile(outputPath, sprintf('of-flow2-graythresh-detail1-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(4));
export_fig(fullfile(outputPath, sprintf('of-flow2-graythresh-detail2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(5));
export_fig(fullfile(outputPath, sprintf('of-flow2-graythresh-detail2-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(6));
close all;

% Set regularisation parameters.
alpha = 0.01;
beta = 0.001;
gamma = 0.005;

% Create folder string with parameters.
folderstr = sprintf('alpha-%.4g-beta-%.4g-gamma-%.4g', alpha, beta, gamma);

% Solve linear system for mass conservation.
timerVal = tic;
c = (Acm + alpha * D + beta * E + gamma * G) \ bcm;
relres = norm((Acm + alpha * D + beta * E + gamma * G) * c - bcm) / norm(bcm);
elapsed = toc(timerVal);
fprintf('Relative residual %e.\n', relres);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Compute flow.
u = bsxfun(@times, full((bfc1')*c), d1) + bsxfun(@times, full((bfc2')*c), d2);

% Plot flow.
plotflow(F, S, ICS, fd1, fd2, u, cmap);

% Save figures.
export_fig(fullfile(outputPath, sprintf('cm-flow3-graythresh-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(1));
export_fig(fullfile(outputPath, sprintf('cm-flow2-graythresh-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(2));
export_fig(fullfile(outputPath, sprintf('cm-flow2-graythresh-detail1-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(3));
export_fig(fullfile(outputPath, sprintf('cm-flow2-graythresh-detail1-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(4));
export_fig(fullfile(outputPath, sprintf('cm-flow2-graythresh-detail2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(5));
export_fig(fullfile(outputPath, sprintf('cm-flow2-graythresh-detail2-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(6));
close all;

%% Use constant one as segmentation.

% Use constant one.
segfh = @(x) ones(size(x, 1), 1);

fprintf('Computing optimaly conditions.\n');
timerVal = tic;
[~, Aof, Acm, D, E, G, bof, bcm] = optcondofcm(Ns, cs{t}, cs{t+1}, X, k, h, xi, weights, gradfx, dtfx, fx1, segfh(fx1), mem);
elapsed = toc(timerVal);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Set regularisation parameters.
alpha = 0.1;
beta = 0;

% Create folder string with parameters.
folderstr = sprintf('alpha-%.4g-beta-%.4g', alpha, beta);

% Solve linear system for optical flow.
timerVal = tic;
c = (Aof + alpha * D + beta * E) \ bof;
relres = norm((Aof + alpha * D + beta * E) * c - bof) / norm(bof);
elapsed = toc(timerVal);
fprintf('Relative residual %e.\n', relres);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Compute flow.
w = bsxfun(@times, full((bfc1')*c), d1) + bsxfun(@times, full((bfc2')*c), d2);

% Plot flow.
plotflow(F, S, ICS, fd1, fd2, w, cmap);

% Save figures.
export_fig(fullfile(outputPath, sprintf('of-flow3-one-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(1));
export_fig(fullfile(outputPath, sprintf('of-flow2-one-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(2));
export_fig(fullfile(outputPath, sprintf('of-flow2-one-detail1-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(3));
export_fig(fullfile(outputPath, sprintf('of-flow2-one-detail1-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(4));
export_fig(fullfile(outputPath, sprintf('of-flow2-one-detail2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(5));
export_fig(fullfile(outputPath, sprintf('of-flow2-one-detail2-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(6));
close all;

% Set regularisation parameters.
alpha = 0.01;
beta = 0.001;
gamma = 0.005;

% Create folder string with parameters.
folderstr = sprintf('alpha-%.4g-beta-%.4g-gamma-%.4g', alpha, beta, gamma);

% Solve linear system for mass conservation.
timerVal = tic;
c = (Acm + alpha * D + beta * E + gamma * G) \ bcm;
relres = norm((Acm + alpha * D + beta * E + gamma * G) * c - bcm) / norm(bcm);
elapsed = toc(timerVal);
fprintf('Relative residual %e.\n', relres);
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Compute flow.
u = bsxfun(@times, full((bfc1')*c), d1) + bsxfun(@times, full((bfc2')*c), d2);

% Plot flow.
plotflow(F, S, ICS, fd1, fd2, u, cmap);

% Save figures.
export_fig(fullfile(outputPath, sprintf('cm-flow3-one-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(1));
export_fig(fullfile(outputPath, sprintf('cm-flow2-one-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(2));
export_fig(fullfile(outputPath, sprintf('cm-flow2-one-detail1-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(3));
export_fig(fullfile(outputPath, sprintf('cm-flow2-one-detail1-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(4));
export_fig(fullfile(outputPath, sprintf('cm-flow2-one-detail2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1', figure(5));
export_fig(fullfile(outputPath, sprintf('cm-flow2-one-detail2-%s-%s-%.3i.png', name, folderstr, frames(t)+1)), '-png', quality, '-transparent', '-a1', figure(6));
close all;

%% Streamlines.

% Set seed points for streamlines.
[xs, ys] = meshgrid(-400:15:400, -400:15:400);

% Set parameters.
stepsize = 1;
maxit = 50;
lineWidth = 1;

% Plot streamlines for flow.
figure(1);
colormap('summer');
hold on;
view(2);
adjust3dplot;
streamlines2(ICS(:, 1:2), w(:, 1:2), [xs(:), ys(:)], stepsize, maxit, 'summer', lineWidth);

%% Basis functions.

% Evaluate basis functions at vertices.
[bfc1, bfc2] = vbasiscompmem(k, h, X, IC, mem);

% Compute support of basis functions.
supp = bfc1 > 0;
supp = full(sum(supp, 1))';

% Plot sum of basis functions.
figure(1);
daspect([1, 1, 1]);
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', supp);
set(gca, 'ZLim', [0, 1]);
set(gca, 'XLim', [-1, 1]);
set(gca, 'YLim', [-1, 1]);
set(gca, 'XTick', -1:0.5:1);
set(gca, 'YTick', -1:0.5:1);
set(gca, 'ZTick', -1:0.5:1);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);
colorbar;
view(3);
