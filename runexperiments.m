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
%frames = 112:115;
frames = 110:120;

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
deg = 600;

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
s = cellfun(@(x) double(im2bw(x, graythresh(x))), fx, 'UniformOutput', false);
s = fx;
s = cellfun(@(x) ones(size(x, 1), 1), fx, 'UniformOutput', false);

% Run through all pair of frames.
for t=1:length(frames)-1
    fprintf('Computing velocity field %i/%i.\n', t, length(frames)-1);
    
    % Compute optimality conditions.
    [~, A, D, E, G, b] = optcondcm(Ns, cs{t}, cs{t+1}, X, k, h, xi, w, gradfx{t}, dtfx{t}, fx{t}, s{t}, mem);

    % Solve linear system.
    tic;
    [ofc{t}, L{t}] = solvesystem(A + alpha * D + beta * E + gamma * G, b, 1e-6, 2000);
    fprintf('GMRES terminated at iteration %i with relative residual %e.\n', L{t}.iter(2), L{t}.relres);
    toc;

end

% Save experiment.
save(fullfile(outputPath, sprintf('%s.mat', name)));