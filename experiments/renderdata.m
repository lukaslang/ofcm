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
timestamp1 = '2017-05-19-14-37-42';

% Set datestring of experiment.
timestamp2 = '2017-05-19-21-38-04';

% Define render quality (set to '-r600' for print quality).
quality = '-r300';

% Load data.
path = fullfile('results', name);
file = fullfile(path, sprintf('%s-data.mat', timestamp1));
load(file);

% Load colormap.
load(fullfile('data', 'cmapblue.mat'));

% Specify max. memory for matrix multiplication.
mem = 8*1024^3;

% Create triangulation for visualisation purpose.
[F, V] = halfsphTriang(7);
[S, rho] = cellfun(@(c) surfsynth(Ns, V, c), cs, 'UniformOutput', false);

% Evaluate data at vertices.
fd = evaldata(f, scale, S, sc, bandwidth, layers);

% Create segmentation.
sfd = cellfun(@(x) segfh(x), fd, 'UniformOutput', false);

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

% Find min and max values of rho.
rhomin = min(min([rho{:}], [], 2));
rhomax = max(max([rho{:}], [], 2));

% Create output folders.
outputPath = fullfile('results', name, timestamp1, timestamp2);
mkdir(outputPath);

% Plot rho.
mkdir(fullfile(outputPath, 'rho2'));
for t=1:length(frames)
    % Plot function rho on the unit sphere.
    figure(1);
    cla;
    daspect([1, 1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), rho{t}, 'EdgeColor', 'none');
    shading interp;
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
    if(rhomax - rhomin > eps)
        caxis([rhomin, rhomax]);
    end
    cbar = findobj(figure(1), 'tag', 'Colorbar');
    set(cbar, 'YTick', 280:20:460);
    set(cbar, 'TickLength', 0.02, 'YColor', [0, 0, 0]);
    view(2);
    export_fig(fullfile(outputPath, 'rho2', sprintf('rho2-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
end
close all;

% Plot data.
mkdir(fullfile(outputPath, 'data2'));
mkdir(fullfile(outputPath, 'data3'));
mkdir(fullfile(outputPath, 'raw3'));
for t=1:length(frames)
    
    % Raw data.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    vol3d('cdata', f{t}, 'XData', scale(1) * [0, size(f{t}, 1)] - sc(1), 'YData', scale(2) * [0, size(f{t}, 2)] - sc(2), 'ZData', scale(3) * [0, size(f{t}, 3)] - sc(3));
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'raw3', sprintf('raw3-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    
    % Raw data and surface with data.
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    export_fig(fullfile(outputPath, 'raw3', sprintf('raw3-%s-surf-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    
    % Add grid.
    surf(Xm{t}, Ym{t}, Zm{t}, ones(nsphere-2, nsphere), 'EdgeColor', [0.7, 0.7, 0.7], 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeAlpha', 1);
    export_fig(fullfile(outputPath, 'raw3', sprintf('raw3-%s-surf-grid-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    
    % Cross-section of narrow band.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    vol3d('cdata', f{t}, 'XData', scale(1) * [0, size(f{t}, 1)] - sc(1), 'YData', scale(2) * [0, size(f{t}, 2)] - sc(2), 'ZData', scale(3) * [0, size(f{t}, 3)] - sc(3));
    view([-90, 0]);
    adjust3dplot;
    set(gca, 'XLim', [0, 50]);
    
    % Define narrow band around surface.
    Sinner = cellfun(@(x) bandwidth(1) * x, S, 'UniformOutput', false);
    Souter = cellfun(@(x) bandwidth(2) * x, S, 'UniformOutput', false);
    trisurf(F, Sinner{t}(:, 1), Sinner{t}(:, 2), Sinner{t}(:, 3), zeros(size(fd{t}, 1), 1), 'EdgeColor', 'green', 'LineWidth', 1.5);
    trisurf(F, Souter{t}(:, 1), Souter{t}(:, 2), Souter{t}(:, 3), zeros(size(fd{t}, 1), 1), 'EdgeColor', 'red', 'LineWidth', 1.5);
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), zeros(size(fd{t}, 1), 1), 'EdgeColor', [0.7, 0.7, 0.7], 'LineWidth', 1.5);
    export_fig(fullfile(outputPath, 'raw3', sprintf('raw3-%s-cross-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    
    % Surface data.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');

    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'data2', sprintf('data2-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    
    % Create a plot with grid.
    view(3);
    surf(Xm{t}, Ym{t}, Zm{t}, ones(nsphere-2, nsphere), 'EdgeColor', [0.7, 0.7, 0.7], 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeAlpha', 1);
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-grid-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    
    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'data2', sprintf('data2-%s-grid-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
end
close all;

% Plot segmentation.
mkdir(fullfile(outputPath, 'seg2'));
for t=1:length(frames)-1
    % Segmentation data.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), sfd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(2);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'seg2', sprintf('seg2-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
end
close all;

% Compute surface normals.
N = cellfun(@(c) surfnormals(Ns, c, xi), cs(1:end-1), 'UniformOutput', false);

% Compute surface velocity.
[~, dtrho] = cellfun(@(x, y) surfsynth(Ns, IC, (y - x) / dt), cs(1:end - 1), cs(2:end), 'UniformOutput', false);
Vs = cellfun(@(x) x .* IC, dtrho, 'UniformOutput', false);
Vsnorm = cellfun(@(x) max(sqrt(sum(x.^2, 2))), Vs, 'UniformOutput', false);
Vsmax = max([Vsnorm{:}]);

% Project surface velocity and compute max norm.
Vsp = cellfun(@(x) projecttoplane(x), Vs, 'UniformOutput', false);
Vspnorm = cellfun(@(x) max(sqrt(sum(x.^2, 2))), Vsp, 'UniformOutput', false);
Vspmax = max([Vspnorm{:}]);

% Compute scalar normal part of surface velocity.
Vsn = cellfun(@(x, y) dot(x, y, 2), Vs, N, 'UniformOutput', false);

% Find min and max values of veln.
Vsnmin = min(min([Vsn{:}], [], 2));
Vsnmax = max(max([Vsn{:}], [], 2));

% Compute signed norm.
signednorm = cellfun(@(x, y) sign(x) .* sqrt(sum(y.^2, 2)), dtrho, Vs, 'UniformOutput', false);

% Find min and max values of signednorm.
signednormmin = min(min([signednorm{:}], [], 2));
signednormmax = max(max([signednorm{:}], [], 2));

% Plot surface velocity.
mkdir(fullfile(outputPath, 'surfvel2'));
mkdir(fullfile(outputPath, 'surfvel3'));
mkdir(fullfile(outputPath, 'surfveln2'));
mkdir(fullfile(outputPath, 'surfveln3'));
mkdir(fullfile(outputPath, 'signednormsurfvel2'));
mkdir(fullfile(outputPath, 'signednormsurfvel3'));
for t=1:length(frames)-1

    % Compute colour of projection.
    col = double(squeeze(computeColour(Vsp{t}(:, 1)/Vspmax, Vsp{t}(:, 2)/Vspmax))) ./ 255;

    % Surface velocity.
    figure(1);
    cla;
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
    C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'surfvel3', sprintf('surfvel3-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    delete(C);
    
    % Top view.
    view(2);
    C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
    export_fig(fullfile(outputPath, 'surfvel2', sprintf('surfvel2-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');

    % Normal part.
    figure(1);
    cla;
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', Vsn{t});
    view(3);
    colorbar;
    if(Vsnmax - Vsnmin > eps)
        caxis([Vsnmin, Vsnmax]);
    end
    adjust3dplot;
    export_fig(fullfile(outputPath, 'surfveln3', sprintf('surfveln3-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    
    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'surfveln2', sprintf('surfveln2-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    colorbar('off');
    
    % Signed norm.
    figure(1);
    cla;
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', signednorm{t});
    view(3);
    colorbar;
    if(signednormmax - signednormmin > eps)
        caxis([signednormmin, signednormmax]);
    end
    adjust3dplot;
    export_fig(fullfile(outputPath, 'signednormsurfvel3', sprintf('signednorm3-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    
    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'signednormsurfvel2', sprintf('signednorm2-%s-%.3i.png', name, frames(t))), '-png', quality, '-transparent', '-a1');
    colorbar('off');
end
close all;