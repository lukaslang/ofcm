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
% three files:
%
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]-data.mat
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]-linsys-of.mat
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]-linsys-cm.mat
%
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
path = fullfile(datapath, 'LSM 16.03.2012');
file = fullfile(path, strcat(name, '.lsm'));

% Select frames.
frames = 112:115;
%frames = 110:120;

% Create start date and time.
startdate = datestr(now, 'yyyy-mm-dd-HH-MM-SS');

% Define and create output folder.
outputPath = fullfile('results', name);
mkdir(outputPath);

% Specify max. memory for matrix multiplication.
mem = 5*1024^3;

% Define Gaussian filter.
sigma = 1.5;
hsize = 30;

% Define thresholding filter.
threshold = 80;
area = 100;

% Define surface fitting parameters.
Ns = 0:30;
beta0 = 1e-4;
beta1 = 0.5;
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
%seg = cellfun(@(x) double(im2bw(x, graythresh(x))), fx, 'UniformOutput', false);
seg = fx;
%seg = cellfun(@(x) ones(size(x, 1), 1), fx, 'UniformOutput', false);

% Save data and parameters.
save(fullfile(outputPath, sprintf('%s-data.mat', startdate)), 'name', 'file', 'frames', 'outputPath', 'Ns', 'cs', 'f', 'scale', 'sc', 'bandwidth', 'layers', 'k', 'h', 'X', 'mem', 'beta0', 'beta1', 's', 'dt', 'ref', 'deg', 'sigma', 'hsize', 'threshold', 'seg', '-v7.3');

% Initialise cell arrays for optical flow.
A = cell(length(frames)-1, 1);
D = cell(length(frames)-1, 1);
E = cell(length(frames)-1, 1);
b = cell(length(frames)-1, 1);

% Run through all pair of frames.
for t=1:length(frames)-1
    fprintf('Computing optimality conditions for OF %i/%i.\n', t, length(frames)-1);
    
    % Compute optimality conditions.
    timerVal = tic;
    [~, A{t}, D{t}, E{t}, b{t}] = optcond(Ns, cs{t}, X, k, h, xi, w, gradfx{t}, dtfx{t}, seg{t}, mem);
    elapsed = toc(timerVal);
    fprintf('Elapsed time is %.6f seconds.\n', elapsed);
end

% Save linear systems.
save(fullfile(outputPath, sprintf('%s-linsys-of.mat', startdate)), 'A', 'D', 'E', 'b', '-v7.3');

% Initialise cell arrays for mass conservation.
A = cell(length(frames)-1, 1);
D = cell(length(frames)-1, 1);
E = cell(length(frames)-1, 1);
G = cell(length(frames)-1, 1);
b = cell(length(frames)-1, 1);

% Run through all pair of frames.
for t=1:length(frames)-1
    fprintf('Computing optimality conditions for CM %i/%i.\n', t, length(frames)-1);
    
    % Compute optimality conditions.
    timerVal = tic;
    [~, A{t}, D{t}, E{t}, G{t}, b{t}] = optcondcm(Ns, cs{t}, cs{t+1}, X, k, h, xi, w, gradfx{t}, dtfx{t}, fx{t}, seg{t}, mem);
    elapsed = toc(timerVal);
    fprintf('Elapsed time is %.6f seconds.\n', elapsed);
end

% Save linear systems.
save(fullfile(outputPath, sprintf('%s-linsys-cm.mat', startdate)), 'A', 'D', 'E', 'G', 'b', '-v7.3');