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

% Create segmentation.
%sfd = cellfun(@(x) double(im2bw(x, graythresh(x))), fd, 'UniformOutput', false);
sfd = fd;
%sfd = cellfun(@(x) ones(size(x, 1), 1), fd, 'UniformOutput', false);

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

% Select experiment.
p = 1;

% Extract components.
%for t=1:length(frames)-1
%    % Project and scale.
%    vp{t, p} = projecttoplane(v{t, p});
%end
%vp = cat(3, vp{:, p});
vp = cat(3, v{:, p});
vx = squeeze(vp(:, 1, :));
vy = squeeze(vp(:, 2, :));

% Subtract mean.
vx = vx - mean(vx, 2);
vy = vy - mean(vy, 2);

n = 2^nextpow2(size(vx, 2));

% Compute Fourier transform.
yx = fft(vx', n);
yy = fft(vy', n);

% Plot power.
n = size(yx, 1);
powerx = abs(yx(1:floor(n/2), :)).^2;
powery = abs(yy(1:floor(n/2), :)).^2;
maxfreq = 1/2;
freq = (1:n/2)'/(n/2)*maxfreq;
%figure(1);
%hold on;
%plot(freq, mean(powerx, 2));
%plot(freq, mean(powery, 2));

% Low-pass filter.
mask = zeros(n, 1);
f = 0:n/2;
sigmaf = 1;
cfreq = 2;
amp = 1;
mask(1:n/2+1) = amp*exp(-((f-cfreq)/(2*sigmaf)).^2);
mask(n:-1:n/2+2) = mask(2:n/2);
figure(2);
plot(mask); title('Mask');

% Multiply with the mask
yxm = yx .* mask;
yym = yy .* mask;

% Plot power of mask.
n = size(yx, 1);
powerx = abs(yxm(1:floor(n/2), :)).^2;
powery = abs(yym(1:floor(n/2), :)).^2;
maxfreq = 1/2;
freq = (1:n/2)'/(n/2)*maxfreq;
%figure(3);
%hold on;
%plot(freq, mean(powerx, 2));
%plot(freq, mean(powery, 2));

% Apply inverse Fourier transform.
vxi = ifft(yxm)';
vyi = ifft(yym)';

len = cell(length(frames)-1, 1);
nmax = -inf;
for t=1:length(frames)-1
    len{t} = sqrt(vxi(:, t).^2 + vyi(:, t).^2);
    % Compute colour space scaling.
    nmax = max(nmax, max(sqrt(len{t})));
end

% Plot filtered velocity fields.
for t=1:length(frames)-1
    % Compute colour of projection.
    col = double(squeeze(computeColour(vxi(:, t)/nmax, vyi(:, t)/nmax))) ./ 255;
    
    disp(t);
    figure(4);
    subplot(1, 2, 1);
    cla;
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
    view(2);
    adjust3dplot;
    
    subplot(1, 2, 2);
    cla;
    colormap(cmap);
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);
    adjust3dplot;
    
    drawnow();
end