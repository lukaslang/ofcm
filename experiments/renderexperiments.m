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
timestamp1 = '2017-05-16-13-50-15';

% Set datestring of experiment.
timestamp2 = '2017-05-16-14-35-42';

% Load data.
path = fullfile('results', name);
file = fullfile(path, sprintf('%s-data.mat', timestamp1));
load(file);

% Load colormap.
load(fullfile('data', 'cmapblue.mat'));

% Specify max. memory for matrix multiplication.
mem = 3*1024^3;

% Figure parameters.
zlim = [0, 400];
xlim = [-400, 400];
ylim = [-400, 400];
ztick = -450:150:450;
xtick = -450:150:450;
ytick = -450:150:450;


% Create triangulation for visualisation purpose.
[F, V] = halfsphTriang(7);
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
export_fig(fullfile(outputPath, 'colourwheel.png'), '-png', '-r600', '-a1', '-transparent');
close all;

% Plot rho.
mkdir(fullfile(outputPath, 'rho2'));
mkdir(fullfile(outputPath, 'rho3'));
for t=1:length(frames)-1
    % Plot function rho on the unit sphere.
    figure(1);
    cla;
    axis square;
    daspect([1, 1, 1]);
    hold on;
    if(rhomax - rhomin > eps)
        caxis([rhomin rhomax]);
    end
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
    export_fig(fullfile(outputPath, 'rho3', sprintf('rho3-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');

    % Rotate by pi.
    [az, el] = view;
    view(az + 180, el);
    export_fig(fullfile(outputPath, 'rho3', sprintf('rho3-%s-rotated-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');

    % Top view.
    colorbar;
    cbar = findobj(figure(1), 'tag', 'Colorbar');
    set(cbar, 'YTick', 280:20:460);
    set(cbar, 'TickLength', 0.02, 'YColor', [0, 0, 0]);
    view(2);
    export_fig(fullfile(outputPath, 'rho2', sprintf('rho2-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
end
close all;

% Plot data.
mkdir(fullfile(outputPath, 'data2'));
mkdir(fullfile(outputPath, 'data3'));
mkdir(fullfile(outputPath, 'raw3'));
for t=1:length(frames)-1
    
    % Raw data.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    vol3d('cdata', f{t}, 'XData', scale(2) * [0, 512] - sc(1), 'YData', scale(1) * [0, 512] - sc(2), 'ZData', scale(3) * [0, 44] - sc(3));
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'raw3', sprintf('raw3-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
    
    % Raw data and surface with data.
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    export_fig(fullfile(outputPath, 'raw3', sprintf('raw3-%s-surf-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
    
    % Surface data.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');

    % Rotate by pi.
    [az, el] = view;
    view(az + 180, el);
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-rotated-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');

    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'data2', sprintf('data2-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
    
    % Create a plot with grid.
    view(3);
    surf(Xm{t}, Ym{t}, Zm{t}, ones(nsphere-2, nsphere), 'EdgeColor', [0.7, 0.7, 0.7], 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeAlpha', 1);
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-grid-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
    
    % Rotate by pi.
    [az, el] = view;
    view(az + 180, el);
    export_fig(fullfile(outputPath, 'data3', sprintf('data3-%s-grid-rotated-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');

    % Top view.
    view(2);
    export_fig(fullfile(outputPath, 'data2', sprintf('data2-%s-grid-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
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
    export_fig(fullfile(outputPath, 'seg2', sprintf('seg2-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
end
close all;

% Compute surface normals.
N = cellfun(@(c) surfnormals(Ns, c, xi), cs(1:end-1), 'UniformOutput', false);

% Compute surface velocity.
[~, dtrho] = cellfun(@(x, y) surfsynth(Ns, IC, (y - x) / dt), cs(1:end - 1), cs(2:end), 'UniformOutput', false);
vel = cellfun(@(x) x .* IC, dtrho, 'UniformOutput', false);

% Compute scalar normal part of surface velocity.
veln = cellfun(@(x, y) dot(x, y, 2), vel, N, 'UniformOutput', false);

% Compute signed norm.
signednorm = cellfun(@(x, y) sign(x) .* sqrt(sum(y.^2, 2)), dtrho, vel, 'UniformOutput', false);

% Open a file to save surface velocity.
fid = fopen(fullfile(outputPath, 'surfvelocities.txt'), 'w');
    
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
    fprintf(fid, 'Frame %i: V=%.4f, Vn=%.4f, proj(V)=%.4f\n', t, max(sqrt(sum(vel{t}.^2, 2))), max(sqrt(sum(veln{t}.^2, 2))), nmax);
    
    % Compute colour of projection.
    col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

    % Surface velocity.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
    C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'surfvel3', sprintf('surfvel3-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
    delete(C);
    
    % Top view.
    view(2);
    C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
    export_fig(fullfile(outputPath, 'surfvel2', sprintf('surfvel2-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');

    % Normal part.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', veln{t});
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'surfveln3', sprintf('surfveln3-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
    
    % Top view.
    view(2);
    colorbar;
    export_fig(fullfile(outputPath, 'surfveln2', sprintf('surfveln2-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
    
    % Signed norm.
    figure(1);
    cla;
    colormap(cmap);
    hold on;
    trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', signednorm{t});
    view(3);
    adjust3dplot;
    export_fig(fullfile(outputPath, 'signednormsurfvel3', sprintf('signednorm3-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
    
    % Top view.
    view(2);
    colorbar;
    export_fig(fullfile(outputPath, 'signednormsurfvel2', sprintf('signednorm2-%s-%i-600dpi.png', name, t)), '-png', '-r600', '-transparent', '-a1');
end
close all;

% Close file.
fclose(fid);

% Open a file to save flow scaling.
fid = fopen(fullfile(outputPath, 'velocities-of.txt'), 'w');

% Run through all pairs of frames.
mkdir(fullfile(outputPath, 'vec2-of'));
mkdir(fullfile(outputPath, 'vec3-of'));
mkdir(fullfile(outputPath, 'flow2-of'));
mkdir(fullfile(outputPath, 'flow3-of'));
mkdir(fullfile(outputPath, 'motion2-of'));
mkdir(fullfile(outputPath, 'motion3-of'));
for t=1:length(frames)-1
    fprintf('Rendering frame %i/%i.\n', t, length(frames)-1);

     % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-of-%.3i.mat', timestamp2, t));
    load(file);
    
    % Run through all parameter configurations.
    for p=1:length(c)
        fprintf('Parameter setting %i/%i.\n', p, length(c));
        
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
        export_fig(fullfile(outputPath, 'vec3-of', sprintf('vec3-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');

        % Top view.
        view(2);
        export_fig(fullfile(outputPath, 'vec2-of', sprintf('vec2-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        
        % Project and scale flow.
        up = projecttoplane(v);

        % Compute colour space scaling.
        nmax = max(sqrt(sum(up.^2, 2)));

        % Save flow scaling.
        fprintf(fid, 'Frame %i: w=%.4f, proj(w)=%.4, ', t, max(sqrt(sum(v.^2, 2))), nmax);
    
        % Compute colour of projection.
        col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

        % Optical flow.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'flow3-of', sprintf('flow3-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        delete(C);

        % Rotate by pi.
        [az, el] = view;
        view(az + 180, el);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'flow3-of', sprintf('flow3-%s-setting-%i-rotated-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'flow2-of', sprintf('flow2-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        
        % Recover total velocity.
        up = projecttoplane(vel{t} + v);

        % Compute colour space scaling.
        nmax = max(sqrt(sum(up.^2, 2)));

        % Save flow scaling.
        fprintf(fid, 'U=%.4f, proj(U)=%.4\n', t, max(sqrt(sum((vel{t} + v).^2, 2))), nmax);
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

        % Total velocity.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'motion3-of', sprintf('motion3-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        delete(C);

        % Rotate by pi.
        [az, el] = view;
        view(az + 180, el);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'motion3-of', sprintf('motion3-%s-setting-%i-rotated-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'motion2-of', sprintf('motion2-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
    end
end
close all;

% Close file.
fclose(fid);

% Open a file to save flow scaling.
fid = fopen(fullfile(outputPath, 'velocities-cm.txt'), 'w');

% Run through all pairs of frames.
mkdir(fullfile(outputPath, 'vec2-cm'));
mkdir(fullfile(outputPath, 'vec3-cm'));
mkdir(fullfile(outputPath, 'flow2-cm'));
mkdir(fullfile(outputPath, 'flow3-cm'));
mkdir(fullfile(outputPath, 'motion2-cm'));
mkdir(fullfile(outputPath, 'motion3-cm'));
for t=1:length(frames)-1
    fprintf('Rendering frame %i/%i.\n', t, length(frames)-1);

     % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-cm-%.3i.mat', timestamp2, t));
    load(file);

    % Run through all parameter configurations.
    for p=1:length(c)
        fprintf('Parameter setting %i/%i.\n', p, length(c));
        
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
        export_fig(fullfile(outputPath, 'vec3-cm', sprintf('vec3-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');

        % Top view.
        view(2);
        export_fig(fullfile(outputPath, 'vec2-cm', sprintf('vec2-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        
        % Project and scale flow.
        up = projecttoplane(v);

        % Compute colour space scaling.
        nmax = max(sqrt(sum(up.^2, 2)));

        % Save flow scaling.
        fprintf(fid, 'Frame %i: u=%.4f, proj(u)=%.4, ', t, max(sqrt(sum(v.^2, 2))), nmax);
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

        % Optical flow.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'flow3-cm', sprintf('flow3-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        delete(C);

        % Rotate by pi.
        [az, el] = view;
        view(az + 180, el);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'flow3-cm', sprintf('flow3-%s-setting-%i-rotated-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'flow2-cm', sprintf('flow2-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        
        % Recover total velocity.
        up = projecttoplane(veln{t} .* N{t} + v);

        % Compute colour space scaling.
        nmax = max(sqrt(sum(up.^2, 2)));

        % Save flow scaling.
        fprintf(fid, 'U=%.4f, proj(U)=%.4\n', t, max(sqrt(sum((veln{t} .* N{t} + v).^2, 2))), nmax);
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

        % Total velocity.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'motion3-cm', sprintf('motion3-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        delete(C);

        % Rotate by pi.
        [az, el] = view;
        view(az + 180, el);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'motion3-cm', sprintf('motion3-%s-setting-%i-rotated-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(250:449, -449:-250, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'motion2-cm', sprintf('motion2-%s-setting-%i-%i-600dpi.png', name, p, t)), '-png', '-r600', '-transparent', '-a1');
    end
end
close all;

% Close file.
fclose(fid);

function adjust3dplot
    axis square;
    daspect([1, 1, 1]);
    set(gca, 'ZLim', zlim);
    set(gca, 'XLim', xlim);
    set(gca, 'YLim', ylim);
    set(gca, 'XTick', xtick);
    set(gca, 'YTick', ytick);
    set(gca, 'ZTick', ztick);
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 14);
end