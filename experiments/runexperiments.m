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
% the linear systems, and outputs files:
%
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]/[yyyy-mm-dd-HH-MM-SS]-coeff-of-000.mat
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]/[yyyy-mm-dd-HH-MM-SS]-coeff-cm-000.mat
%
% where the first datestring is from the files prepareexperiments.m:
%
%   results/[name]/[yyyy-mm-dd-HH-MM-SS]-linsys-000.mat
%
% and the second datestring is generated.
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
timestamp = '2017-05-11-23-25-02';

% Define folder.
path = fullfile('results', name);

% Create output folder.
mkdir(fullfile(path, timestamp));

% Create start date and time.
startdate = datestr(now, 'yyyy-mm-dd-HH-MM-SS');

% Set solver parameters.
tolSolver = 1e-6;
iterSolver = 2000;

% Load metadata.
load(fullfile(path, sprintf('%s-data.mat', timestamp)), 'frames');

%% Run experiments for optical flow.

% Set regularisation parameters.
alpha = {0.001, 0.01, 0.1, 1};
beta = {0.001, 0.001, 0.001, 0.001};

% Initialise arrays.
c = cell(length(alpha), 1);
L = cell(length(alpha), 1);

% Run through all pair of frames.
for t=1:length(frames)-1
    fprintf('Solving linear system for OF %i/%i.\n', t, length(frames)-1);
    
    % Load linear systems.
    load(fullfile(path, sprintf('%s-linsys-%.3i.mat', timestamp, t)), 'Aof', 'D', 'E', 'bof');
    
    % Run through parameter settings.
    for p=1:length(alpha)
        % Solve linear system.
        timerVal = tic;
        [c{p}, L{p}] = solvesystem(Aof + alpha{p} * D + beta{p} * E, bof, tolSolver, iterSolver);
        fprintf('GMRES terminated at iteration %i with relative residual %e.\n', L{t}.iter(2), L{t}.relres);
        elapsed = toc(timerVal);
        fprintf('Elapsed time is %.6f seconds.\n', elapsed);
    end
    
    % Save experiments.
    save(fullfile(path, timestamp, sprintf('%s-coeff-of-%.3i.mat', startdate, t)), 'c', 'L', 'alpha', 'beta', '-v7.3');
end

%% Run experiments for mass conservation.

% Set regularisation parameters.
alpha = {0.001, 0.01, 0.1, 1};
beta = {0.001, 0.001, 0.001, 0.001};
gamma = {0.1, 0.1, 0.1, 0.1};

% Initialise arrays.
c = cell(length(alpha), 1);
L = cell(length(alpha), 1);

% Run through all pair of frames.
for t=1:length(frames)-1
    fprintf('Solving linear system for CM %i/%i.\n', t, length(frames)-1);
    
    % Load linear systems.
    load(fullfile(path, sprintf('%s-linsys-%.3i.mat', timestamp, t)), 'Acm', 'D', 'E', 'G', 'bcm');
    
    % Run through parameter settings.
    for p=1:length(alpha)
        % Solve linear system.
        timerVal = tic;
        [c{p}, L{p}] = solvesystem(Acm + alpha{p} * D + beta{p} * E + gamma{p} * G, bcm, tolSolver, iterSolver);
        fprintf('GMRES terminated at iteration %i with relative residual %e.\n', L{t}.iter(2), L{t}.relres);
        elapsed = toc(timerVal);
        fprintf('Elapsed time is %.6f seconds.\n', elapsed);
    end
    
    % Save experiments.
    save(fullfile(path, timestamp, sprintf('%s-coeff-cm-%.3i.mat', startdate, t)), 'c', 'L', 'alpha', 'beta', '-v7.3');
end