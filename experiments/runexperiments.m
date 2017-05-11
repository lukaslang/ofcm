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

% This script loads optimality conditions, computes coefficients by solving
% the linear systems, and outputs two files:
%
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]/[yyyy-mm-dd-HH-MM-SS]-coeff-of.mat
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]/[yyyy-mm-dd-HH-MM-SS]-coeff-cm.mat
%
% where the first datestring is from the files prepareexperiments.m:
%
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]-linsys-of.mat
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]-linsys-cm.mat
%
% and the second datestring is generated.
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
timestamp = '2017-05-11-14-52-57';

% Define folder.
path = fullfile('results', name);

% Create output folder.
mkdir(fullfile(path, timestamp));

% Create start date and time.
startdate = datestr(now, 'yyyy-mm-dd-HH-MM-SS');

% Set solver parameters.
tolSolver = 1e-6;
iterSolver = 2000;

% Load linear systems.
load(fullfile(path, sprintf('%s-linsys.mat', timestamp)));

%% Run experiments for optical flow.

% Set regularisation parameters.
alpha = {0.001, 0.01, 0.1, 1};
beta = {0.001, 0.001, 0.001, 0.001};

% Initialise arrays.
c = cell(size(Aof, 1), length(alpha));
L = cell(size(Aof, 1), length(alpha));

% Run through all pair of frames.
for t=1:size(Aof, 1)
    fprintf('Solving linear system for OF %i/%i.\n', t, size(Aof, 1));
    
    for p=1:length(alpha)
        % Solve linear system.
        timerVal = tic;
        [c{t, p}, L{t, p}] = solvesystem(Aof{t} + alpha{p} * D{t} + beta{p} * E{t}, bof{t}, tolSolver, iterSolver);
        fprintf('GMRES terminated at iteration %i with relative residual %e.\n', L{t}.iter(2), L{t}.relres);
        elapsed = toc(timerVal);
        fprintf('Elapsed time is %.6f seconds.\n', elapsed);
    end
end

% Save experiments.
save(fullfile(path, timestamp, sprintf('%s-coeff-of.mat', startdate)), 'c', 'L', 'alpha', 'beta', '-v7.3');

%% Run experiments for mass conservation.

% Set regularisation parameters.
alpha = {0.001, 0.01, 0.1, 1};
beta = {0.001, 0.001, 0.001, 0.001};
gamma = {0.1, 0.1, 0.1, 0.1};

% Initialise arrays.
c = cell(size(Acm, 1), length(alpha));
L = cell(size(Acm, 1), length(alpha));

% Run through all pair of frames.
for t=1:size(Acm, 1)
    fprintf('Solving linear system for CM %i/%i.\n', t, size(Acm, 1));
    
    for p=1:length(alpha)
        % Solve linear system.
        timerVal = tic;
        [c{t, p}, L{t, p}] = solvesystem(Acm{t} + alpha{p} * D{t} + beta{p} * E{t} + gamma{p} * G{t}, bcm{t}, tolSolver, iterSolver);
        fprintf('GMRES terminated at iteration %i with relative residual %e.\n', L{t}.iter(2), L{t}.relres);
        elapsed = toc(timerVal);
        fprintf('Elapsed time is %.6f seconds.\n', elapsed);
    end
end

% Save experiments.
save(fullfile(path, timestamp, sprintf('%s-coeff-cm.mat', startdate)), 'c', 'L', 'alpha', 'beta', '-v7.3');