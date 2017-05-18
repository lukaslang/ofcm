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
timestamp1 = '2017-05-18-14-19-31';

% Set datestring of experiment.
timestamp2 = '2017-05-18-14-37-56';

% Define render quality.
quality = '-r100';

% Load data.
path = fullfile('results', name);
file = fullfile(path, sprintf('%s-data.mat', timestamp1));
load(file);

% Load colormap.
load(fullfile('data', 'cmapblue.mat'));

% Specify max. memory for matrix multiplication.
mem = 3*1024^3;

% Create triangulation for visualisation purpose.
[F, V] = halfsphTriang(6);
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
rhomin = min(min([rho{:}]));
rhomax = max(max([rho{:}]));

% Create output folders.
outputPath = fullfile('results', name, timestamp1, timestamp2);
mkdir(outputPath);

% Plot 2D colourwheel.
figure(1);
cw = colourwheelbg;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(2);
export_fig(fullfile(outputPath, 'colourwheel.png'), '-png', quality, '-a1', '-transparent');
close all;

% Plot rho.
mkdir(fullfile(outputPath, 'rho2'));
mkdir(fullfile(outputPath, 'rho3'));
for t=1:length(frames)
    % Plot function rho on the unit sphere.
    figure(1);
    cla;
    daspect([1, 1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), rho{t}, 'EdgeColor', 'none');
    shading interp;
    view(3);
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
    export_fig(fullfile(outputPath, 'rho3', sprintf('rho3-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');

    % Rotate by pi.
    [az, el] = view;
    view(az + 180, el);
    export_fig(fullfile(outputPath, 'rho3', sprintf('rho3-%s-rotated-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');

    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'rho2', sprintf('rho2-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
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
    export_fig(fullfile(outputPath, 'raw3', sprintf('raw3-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
    
    % Raw data and surface with data.
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    export_fig(fullfile(outputPath, 'raw3', sprintf('raw3-%s-surf-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
    
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
    export_fig(fullfile(outputPath, 'raw3', sprintf('raw3-%s-cross-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
    
    % Surface data.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');

    % Rotate by pi.
    [az, el] = view;
    view(az + 180, el);
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-rotated-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');

    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'data2', sprintf('data2-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
    
    % Create a plot with grid.
    view(3);
    surf(Xm{t}, Ym{t}, Zm{t}, ones(nsphere-2, nsphere), 'EdgeColor', [0.7, 0.7, 0.7], 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeAlpha', 1);
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-grid-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
    
    % Rotate by pi.
    [az, el] = view;
    view(az + 180, el);
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-grid-rotated-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');

    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'data2', sprintf('data2-%s-grid-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
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
    export_fig(fullfile(outputPath, 'seg2', sprintf('seg2-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
end
close all;

% Compute surface normals.
N = cellfun(@(c) surfnormals(Ns, c, xi), cs(1:end-1), 'UniformOutput', false);

% Compute surface velocity.
[~, dtrho] = cellfun(@(x, y) surfsynth(Ns, IC, (y - x) / dt), cs(1:end - 1), cs(2:end), 'UniformOutput', false);
vel = cellfun(@(x) x .* IC, dtrho, 'UniformOutput', false);

% Compute scalar normal part of surface velocity.
veln = cellfun(@(x, y) dot(x, y, 2), vel, N, 'UniformOutput', false);

% Find min and max values of veln.
velnmin = min(min([veln{:}]));
velnmax = max(max([veln{:}]));

% Compute signed norm.
signednorm = cellfun(@(x, y) sign(x) .* sqrt(sum(y.^2, 2)), dtrho, vel, 'UniformOutput', false);

% Find min and max values of signednorm.
signednormmin = min(min([signednorm{:}]));
signednormmax = max(max([signednorm{:}]));

% Open a file to save surface velocity.
fid = fopen(fullfile(outputPath, 'velocities-surface.txt'), 'w');
    
% Plot surface velocity.
mkdir(fullfile(outputPath, 'surfvel2'));
mkdir(fullfile(outputPath, 'surfvel3'));
mkdir(fullfile(outputPath, 'surfveln2'));
mkdir(fullfile(outputPath, 'surfveln3'));
mkdir(fullfile(outputPath, 'signednormsurfvel2'));
mkdir(fullfile(outputPath, 'signednormsurfvel3'));
for t=1:length(frames)-1
    
    % Project and scale.
    up = projecttoplane(vel{t});

    % Compute colour space scaling.
    nmax = max(sqrt(sum(up.^2, 2)));

    % Save flow scaling.
    fprintf(fid, 'Frame %.3i: V=%.4g,\t Vn=%.4g,\t proj(V)=%.4g\n', t, max(sqrt(sum(vel{t}.^2, 2))), max(sqrt(sum(veln{t}.^2, 2))), nmax);
    
    % Compute colour of projection.
    col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

    % Surface velocity.
    figure(1);
    cla;
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
    C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'surfvel3', sprintf('surfvel3-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
    delete(C);
    
    % Top view.
    view(2);
    C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
    export_fig(fullfile(outputPath, 'surfvel2', sprintf('surfvel2-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');

    % Normal part.
    figure(1);
    cla;
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', veln{t});
    view(3);
    colorbar;
    if(velnmax - velnmin > eps)
        caxis([velnmin, velnmax]);
    end
    adjust3dplot;
    export_fig(fullfile(outputPath, 'surfveln3', sprintf('surfveln3-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
    
    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'surfveln2', sprintf('surfveln2-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
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
    export_fig(fullfile(outputPath, 'signednormsurfvel3', sprintf('signednorm3-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
    
    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'signednormsurfvel2', sprintf('signednorm2-%s-%.3i.png', name, t)), '-png', quality, '-transparent', '-a1');
    colorbar('off');
end
close all;

% Close file.
fclose(fid);

% Compute surface divergence of basis functions.
bfcdiv = cellfun(@(c) surfdivmem(Ns, c, xi, k, h, X, mem), cs, 'UniformOutput', false);

% Precompute surface divergence of flow.
for t=1:length(frames)-1
    % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-of-%.3i.mat', timestamp2, t));
    load(file);
    
    for p=1:length(c)
        sdiv{t, p} = bfcdiv{t}'*c{p};
    end
end

% Find min and max values of surface divergence for each parameter setting.
for p=1:length(c)
    sdivmin{p} = min(min([sdiv{:, p}]));
    sdivmax{p} = max(max([sdiv{:, p}]));
end

% Open a file to save flow scaling.
fid = fopen(fullfile(outputPath, 'velocities-of.txt'), 'w');

% Run through all pairs of frames.
mkdir(fullfile(outputPath, 'of-vec2'));
mkdir(fullfile(outputPath, 'of-vec3'));
mkdir(fullfile(outputPath, 'of-flow2'));
mkdir(fullfile(outputPath, 'of-flow3'));
mkdir(fullfile(outputPath, 'of-motion2'));
mkdir(fullfile(outputPath, 'of-motion3'));
mkdir(fullfile(outputPath, 'of-div2'));
mkdir(fullfile(outputPath, 'of-div3'));
for t=1:length(frames)-1
    fprintf('Rendering frame %.3i/%.3i.\n', t, length(frames)-1);

     % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-of-%.3i.mat', timestamp2, t));
    load(file);
    
    % Run through all parameter configurations.
    for p=1:length(c)
        fprintf('Parameter setting %.3i/%.3i.\n', p, length(c));
        
        % Create folder.
        folderstr = sprintf('alpha-%.4f-beta-%.4f', alpha{p}, beta{p});
        mkdir(fullfile(outputPath, 'of-vec2', folderstr));
        mkdir(fullfile(outputPath, 'of-vec3', folderstr));
        mkdir(fullfile(outputPath, 'of-flow2', folderstr));
        mkdir(fullfile(outputPath, 'of-flow3', folderstr));
        mkdir(fullfile(outputPath, 'of-motion2', folderstr));
        mkdir(fullfile(outputPath, 'of-motion3', folderstr));
        mkdir(fullfile(outputPath, 'of-div2', folderstr));
        mkdir(fullfile(outputPath, 'of-div3', folderstr));
        
        % Compute pushforward of basis functions.
        v = bsxfun(@times, full((bfc1')*c{p}), d1{t}) + bsxfun(@times, full((bfc2')*c{p}), d2{t});

        % Optical flow vectors.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(3);
        adjust3dplot;
        quiver3(ICS{t}(:, 1), ICS{t}(:, 2), ICS{t}(:, 3), v(:, 1), v(:, 2), v(:, 3), 0, 'r');
        export_fig(fullfile(outputPath, 'of-vec3', folderstr, sprintf('vec3-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Top view.
        view(2);
        export_fig(fullfile(outputPath, 'of-vec2', folderstr, sprintf('vec2-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Optical flow vectors scaled.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(3);
        adjust3dplot;
        quiver3(ICS{t}(:, 1), ICS{t}(:, 2), ICS{t}(:, 3), v(:, 1), v(:, 2), v(:, 3), 1, 'r');
        export_fig(fullfile(outputPath, 'of-vec3', folderstr, sprintf('vec3-%s-%s-scaled-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Top view.
        view(2);
        export_fig(fullfile(outputPath, 'of-vec2', folderstr, sprintf('vec2-%s-%s-scaled-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Project and scale flow.
        up = projecttoplane(v);

        % Compute colour space scaling.
        nmax = max(sqrt(sum(up.^2, 2)));

        % Save flow scaling.
        fprintf(fid, 'Frame %.3i, alpha=%.4f, beta=%.4f: w=%.4g,\t proj(w)=%.4g,\t ', t, alpha{p}, beta{p}, max(sqrt(sum(v.^2, 2))), nmax);
    
        % Compute colour of projection.
        col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

        % Optical flow.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'of-flow3', folderstr, sprintf('flow3-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        delete(C);

        % Rotate by pi.
        [az, el] = view;
        view(az + 180, el);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'of-flow3', folderstr, sprintf('flow3-%s-%s-rotated-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'of-flow2', folderstr, sprintf('flow2-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Recover total velocity.
        up = projecttoplane(vel{t} + v);

        % Compute colour space scaling.
        nmax = max(sqrt(sum(up.^2, 2)));

        % Save flow scaling.
        fprintf(fid, 'V=%.4g,\t V+w=%.4g,\t proj(V+w)=%.4g\n', max(sqrt(sum(vel{t}.^2, 2))), max(sqrt(sum((vel{t} + v).^2, 2))), nmax);
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

        % Total velocity.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'of-motion3', folderstr, sprintf('motion3-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        delete(C);

        % Rotate by pi.
        [az, el] = view;
        view(az + 180, el);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'of-motion3', folderstr, sprintf('motion3-%s-%s-rotated-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'of-motion2', folderstr, sprintf('motion2-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Plot divergence.
        figure(1);
        cla;
        colormap default;
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sdiv{t, p});
        view(3);
        adjust3dplot;
        colorbar;
        if(sdivmax{p} - sdivmin{p} > eps)
            caxis([sdivmin{p}, sdivmax{p}]);
        end
        export_fig(fullfile(outputPath, 'of-div3', folderstr, sprintf('div3-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');

        % Top view.
        view(2);
        export_fig(fullfile(outputPath, 'of-div2', folderstr, sprintf('div2-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        close all;
    end
end
close all;

% Close file.
fclose(fid);

% Open a file to save flow scaling.
fid = fopen(fullfile(outputPath, 'velocities-cm.txt'), 'w');

% Precompute surface divergence of flow.
for t=1:length(frames)-1
    % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-cm-%.3i.mat', timestamp2, t));
    load(file);
    
    for p=1:length(c)
        sdiv{t, p} = bfcdiv{t}'*c{p};
    end
end

% Find min and max values of surface divergence for each parameter setting.
for p=1:length(c)
    sdivmin{p} = min(min([sdiv{:, p}]));
    sdivmax{p} = max(max([sdiv{:, p}]));
end
    
% Run through all pairs of frames.
mkdir(fullfile(outputPath, 'cm-vec2'));
mkdir(fullfile(outputPath, 'cm-vec3'));
mkdir(fullfile(outputPath, 'cm-flow2'));
mkdir(fullfile(outputPath, 'cm-flow3'));
mkdir(fullfile(outputPath, 'cm-motion2'));
mkdir(fullfile(outputPath, 'cm-motion3'));
mkdir(fullfile(outputPath, 'cm-div2'));
mkdir(fullfile(outputPath, 'cm-div3'));
for t=1:length(frames)-1
    fprintf('Rendering frame %.3i/%.3i.\n', t, length(frames)-1);

    % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-cm-%.3i.mat', timestamp2, t));
    load(file);

    % Run through all parameter configurations.
    for p=1:length(c)
        fprintf('Parameter setting %.3i/%.3i.\n', p, length(c));
        
        % Create folder.
        folderstr = sprintf('alpha-%.4f-beta-%.4f-gamma-%.4f', alpha{p}, beta{p}, gamma{p});
        mkdir(fullfile(outputPath, 'cm-vec2', folderstr));
        mkdir(fullfile(outputPath, 'cm-vec3', folderstr));
        mkdir(fullfile(outputPath, 'cm-flow2', folderstr));
        mkdir(fullfile(outputPath, 'cm-flow3', folderstr));
        mkdir(fullfile(outputPath, 'cm-motion2', folderstr));
        mkdir(fullfile(outputPath, 'cm-motion3', folderstr));
        mkdir(fullfile(outputPath, 'cm-div2', folderstr));
        mkdir(fullfile(outputPath, 'cm-div3', folderstr));
        
        % Compute pushforward of basis functions.
        u = bsxfun(@times, full((bfc1')*c{p}), d1{t}) + bsxfun(@times, full((bfc2')*c{p}), d2{t});

        % Optical flow vectors.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(3);
        adjust3dplot;
        quiver3(ICS{t}(:, 1), ICS{t}(:, 2), ICS{t}(:, 3), u(:, 1), u(:, 2), u(:, 3), 0, 'r');
        export_fig(fullfile(outputPath, 'cm-vec3', folderstr, sprintf('vec3-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');

        % Top view.
        view(2);
        export_fig(fullfile(outputPath, 'cm-vec2', folderstr, sprintf('vec2-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Optical flow vectors scaled.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(3);
        adjust3dplot;
        quiver3(ICS{t}(:, 1), ICS{t}(:, 2), ICS{t}(:, 3), u(:, 1), u(:, 2), u(:, 3), 1, 'r');
        export_fig(fullfile(outputPath, 'cm-vec3', folderstr, sprintf('vec3-%s-%s-scaled-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');

        % Top view.
        view(2);
        export_fig(fullfile(outputPath, 'cm-vec2', folderstr, sprintf('vec2-%s-%s-scaled-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Project and scale flow.
        up = projecttoplane(u);

        % Compute colour space scaling.
        nmax = max(sqrt(sum(up.^2, 2)));

        % Save flow scaling.
        fprintf(fid, 'Frame %.3i, alpha=%.4f, beta=%.4f, gamma=%.4f: u=%.4g,\t proj(u)=%.4g,\t ', t, alpha{p}, beta{p}, gamma{p}, max(sqrt(sum(u.^2, 2))), nmax);
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

        % Optical flow.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'cm-flow3', folderstr, sprintf('flow3-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        delete(C);

        % Rotate by pi.
        [az, el] = view;
        view(az + 180, el);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'cm-flow3', folderstr, sprintf('flow3-%s-%s-rotated-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'cm-flow2', folderstr, sprintf('flow2-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Recover total velocity.
        up = projecttoplane(veln{t} .* N{t} + u);

        % Compute colour space scaling.
        nmax = max(sqrt(sum(up.^2, 2)));

        % Save flow scaling.
        fprintf(fid, 'V*N=%.4g,\t V*N+u=%.4g,\t proj(V*N+u)=%.4g\n', max(sqrt(sum((veln{t} .* N{t}).^2, 2))), max(sqrt(sum((veln{t} .* N{t} + u).^2, 2))), nmax);
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

        % Total velocity.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'cm-motion3', folderstr, sprintf('motion3-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        delete(C);

        % Rotate by pi.
        [az, el] = view;
        view(az + 180, el);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'cm-motion3', folderstr, sprintf('motion3-%s-%s-rotated-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(300:399, -399:-300, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'cm-motion2', folderstr, sprintf('motion2-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        
        % Plot divergence.
        figure(1);
        cla;
        colormap default;
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sdiv{t, p});
        view(3);
        adjust3dplot;
        colorbar;
        if(sdivmax{p} - sdivmin{p} > eps)
            caxis([sdivmin{p}, sdivmax{p}]);
        end
        export_fig(fullfile(outputPath, 'cm-div3', folderstr, sprintf('div3-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');

        % Top view.
        view(2);
        export_fig(fullfile(outputPath, 'cm-div2', folderstr, sprintf('div2-%s-%s-%.3i.png', name, folderstr, t)), '-png', quality, '-transparent', '-a1');
        close all;
    end
end
close all;

% Close file.
fclose(fid);