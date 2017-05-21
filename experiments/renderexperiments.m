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

% Define render quality.
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
rhomin = min(min([rho{:}], [], 2));
rhomax = max(max([rho{:}], [], 2));

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

% Compute surface divergence of basis functions.
bfcdiv = cellfun(@(c) surfdivmem(Ns, c, xi, k, h, X, mem), cs, 'UniformOutput', false);

% Precompute.
for t=1:length(frames)-1
    % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-of-%.3i.mat', timestamp2, frames(t)));
    load(file);
    
    for p=1:length(c)
        % Compute divergence.
        sdiv{t, p} = bfcdiv{t}'*c{p};
        
        % Compute flow.
        w{t, p} = bsxfun(@times, full((bfc1')*c{p}), d1{t}) + bsxfun(@times, full((bfc2')*c{p}), d2{t});
        
        % Compute norm of flow.
        wnorm{t, p} = sqrt(sum(w{t, p}.^2, 2));
        
        % Project and scale flow.
        wp{t, p} = projecttoplane(w{t, p});

        % Compute colour space scaling.
        wpnorm{t, p} = sqrt(sum(wp{t, p}.^2, 2));
        
        % Recover total velocity.
        U{t, p} = Vs{t} + w{t, p};
        Unorm{t, p} = sqrt(sum(U{t, p}.^2, 2));
        
        % Project and scale flow.
        Up{t, p} = projecttoplane(U{t, p});
        Upnorm{t, p} = sqrt(sum(Up{t, p}.^2, 2));
    end
end

% Find min and max values.
for p=1:length(c)
    sdivmin{p} = min(min([sdiv{:, p}], [], 2));
    sdivmax{p} = max(max([sdiv{:, p}], [], 2));
    wmin{p} = min(min([wnorm{:, p}], [], 2));
    wmax{p} = max(max([wnorm{:, p}], [], 2));
    wpmax{p} = max(max([wpnorm{:, p}], [], 2));
    Umax{p} = max(max([Unorm{:, p}], [], 2));
    Upmax{p} = max(max([Upnorm{:, p}], [], 2));
end

% Open a file to save flow scaling.
fid = fopen(fullfile(outputPath, 'velocities-of.txt'), 'w');
fprintf(fid, 'max(V)=%.4g, max(proj(V))=%.4g\n', Vsmax, Vspmax);

% Save flow scaling for each experiment.
for p=1:length(c)
    fprintf(fid, 'alpha=%.4g, beta=%.4g: max(w)=%.4g, max(proj(w))=%.4g, max(V+w)=%.4g, max(proj(V+w))=%.4g, divmax(w)=%.4g\n', alpha{p}, beta{p}, wmax{p}, wpmax{p}, Umax{p}, Upmax{p}, sdivmax{p});
end
fclose(fid);

% Run through all pairs of frames.
mkdir(fullfile(outputPath, 'of-vec2'));
mkdir(fullfile(outputPath, 'of-flow2'));
mkdir(fullfile(outputPath, 'of-flow3'));
mkdir(fullfile(outputPath, 'of-flowmagn2'));
mkdir(fullfile(outputPath, 'of-motion2'));
mkdir(fullfile(outputPath, 'of-motion3'));
mkdir(fullfile(outputPath, 'of-div2'));
for t=1:length(frames)-1
    fprintf('Rendering frame %.3i/%.3i.\n', t, length(frames)-1);

     % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-of-%.3i.mat', timestamp2, frames(t)));
    load(file);
    
    % Run through all parameter configurations.
    for p=1:length(c)
        fprintf('Parameter setting %.3i/%.3i.\n', p, length(c));
        
        % Create folder.
        folderstr = sprintf('alpha-%.4g-beta-%.4g', alpha{p}, beta{p});
        mkdir(fullfile(outputPath, 'of-vec2', folderstr));
        mkdir(fullfile(outputPath, 'of-flow2', folderstr));
        mkdir(fullfile(outputPath, 'of-flow3', folderstr));
        mkdir(fullfile(outputPath, 'of-flowmagn2', folderstr));
        mkdir(fullfile(outputPath, 'of-motion2', folderstr));
        mkdir(fullfile(outputPath, 'of-motion3', folderstr));
        mkdir(fullfile(outputPath, 'of-div2', folderstr));
        
        % Optical flow vectors.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        adjust3dplot;
        quiver3(ICS{t}(:, 1), ICS{t}(:, 2), ICS{t}(:, 3), w{t, p}(:, 1), w{t, p}(:, 2), w{t, p}(:, 3), 0, 'r');
        export_fig(fullfile(outputPath, 'of-vec2', folderstr, sprintf('vec2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(wp{t, p}(:, 1)/wpmax{p}, wp{t, p}(:, 2)/wpmax{p}))) ./ 255;

        % Optical flow.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'of-flow3', folderstr, sprintf('flow3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'of-flow2', folderstr, sprintf('flow2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        
        % Plot magnitude of flow.
        figure(1);
        cla;
        colormap default;
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sqrt(sum(w{t, p}.^2, 2)));
        view(2);
        adjust3dplot;
        colorbar;
        if(wmax{p} - wmin{p} > eps)
            caxis([wmin{p}, wmax{p}]);
        end
        export_fig(fullfile(outputPath, 'of-flowmagn2', folderstr, sprintf('flowmagn2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        close all;

        % Compute colour of projection.
        col = double(squeeze(computeColour(Up{t, p}(:, 1)/Upmax{p}, Up{t, p}(:, 2)/Upmax{p}))) ./ 255;
        
        % Total velocity.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'of-motion3', folderstr, sprintf('motion3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'of-motion2', folderstr, sprintf('motion2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        
        % Plot divergence.
        figure(1);
        cla;
        colormap default;
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sdiv{t, p});
        view(2);
        adjust3dplot;
        colorbar;
        if(sdivmax{p} - sdivmin{p} > eps)
            caxis([sdivmin{p}, sdivmax{p}]);
        end
        export_fig(fullfile(outputPath, 'of-div2', folderstr, sprintf('div2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        close all;
    end
end
close all;

% Precompute.
for t=1:length(frames)-1
    % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-cm-%.3i.mat', timestamp2, frames(t)));
    load(file);
    
    for p=1:length(c)
        % Compute divergence.
        sdiv{t, p} = bfcdiv{t}'*c{p};
        
        % Compute flow.
        u{t, p} = bsxfun(@times, full((bfc1')*c{p}), d1{t}) + bsxfun(@times, full((bfc2')*c{p}), d2{t});
        
        % Compute norm of flow.
        unorm{t, p} = sqrt(sum(u{t, p}.^2, 2));
        
        % Project and scale flow.
        up{t, p} = projecttoplane(u{t, p});

        % Compute colour space scaling.
        upnorm{t, p} = sqrt(sum(up{t, p}.^2, 2));
        
        % Recover total velocity.
        U{t, p} = Vsn{t} .* N{t} + u{t, p};
        Unorm{t, p} = sqrt(sum(U{t, p}.^2, 2));
        
        % Project and scale flow.
        Up{t, p} = projecttoplane(U{t, p});
        Upnorm{t, p} = sqrt(sum(Up{t, p}.^2, 2));
    end
end

% Find min and max values.
for p=1:length(c)
    sdivmin{p} = min(min([sdiv{:, p}], [], 2));
    sdivmax{p} = max(max([sdiv{:, p}], [], 2));
    umin{p} = min(min([unorm{:, p}], [], 2));
    umax{p} = max(max([unorm{:, p}], [], 2));
    upmax{p} = max(max([upnorm{:, p}], [], 2));
    Umax{p} = max(max([Unorm{:, p}], [], 2));
    Upmax{p} = max(max([Upnorm{:, p}], [], 2));
end

% Open a file to save flow scaling.
fid = fopen(fullfile(outputPath, 'velocities-cm.txt'), 'w');
fprintf(fid, 'max(V)=%.4g, max(proj(V))=%.4g\n', Vsmax, Vspmax);

% Save flow scaling for each experiment.
for p=1:length(c)
    fprintf(fid, 'alpha=%.4g, beta=%.4g: max(u)=%.4g, max(proj(u))=%.4g, max(VN+u)=%.4g, max(proj(VN+u))=%.4g, divmax(u)=%.4g\n', alpha{p}, beta{p}, umax{p}, upmax{p}, Umax{p}, Upmax{p}, sdivmax{p});
end
fclose(fid);
    
% Run through all pairs of frames.
mkdir(fullfile(outputPath, 'cm-vec2'));
mkdir(fullfile(outputPath, 'cm-flow2'));
mkdir(fullfile(outputPath, 'cm-flow3'));
mkdir(fullfile(outputPath, 'cm-flowmagn2'));
mkdir(fullfile(outputPath, 'cm-motion2'));
mkdir(fullfile(outputPath, 'cm-motion3'));
mkdir(fullfile(outputPath, 'cm-div2'));
for t=1:length(frames)-1
    fprintf('Rendering frame %.3i/%.3i.\n', t, length(frames)-1);

    % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-cm-%.3i.mat', timestamp2, frames(t)));
    load(file);

    % Run through all parameter configurations.
    for p=1:length(c)
        fprintf('Parameter setting %.3i/%.3i.\n', p, length(c));
        
        % Create folder.
        folderstr = sprintf('alpha-%.4g-beta-%.4g-gamma-%.4g', alpha{p}, beta{p}, gamma{p});
        mkdir(fullfile(outputPath, 'cm-vec2', folderstr));
        mkdir(fullfile(outputPath, 'cm-flow2', folderstr));
        mkdir(fullfile(outputPath, 'cm-flow3', folderstr));
        mkdir(fullfile(outputPath, 'cm-flowmagn2', folderstr));
        mkdir(fullfile(outputPath, 'cm-motion2', folderstr));
        mkdir(fullfile(outputPath, 'cm-motion3', folderstr));
        mkdir(fullfile(outputPath, 'cm-div2', folderstr));
        
        % Optical flow vectors.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), fd{t}, 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        adjust3dplot;
        quiver3(ICS{t}(:, 1), ICS{t}(:, 2), ICS{t}(:, 3), u{t, p}(:, 1), u{t, p}(:, 2), u{t, p}(:, 3), 0, 'r');
        export_fig(fullfile(outputPath, 'cm-vec2', folderstr, sprintf('vec2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(up{t, p}(:, 1)/upmax{p}, up{t, p}(:, 2)/upmax{p}))) ./ 255;

        % Optical flow.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'cm-flow3', folderstr, sprintf('flow3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'cm-flow2', folderstr, sprintf('flow2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        
        % Plot magnitude of flow.
        figure(1);
        cla;
        colormap default;
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sqrt(sum(u{t, p}.^2, 2)));
        view(2);
        adjust3dplot;
        colorbar;
        if(umax{p} - umin{p} > eps)
            caxis([umin{p}, umax{p}]);
        end
        export_fig(fullfile(outputPath, 'cm-flowmagn2', folderstr, sprintf('flowmagn2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        close all;
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(Up{t, p}(:, 1)/Upmax{p}, Up{t, p}(:, 2)/Upmax{p}))) ./ 255;

        % Total velocity.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
        C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        view(3);
        adjust3dplot;
        export_fig(fullfile(outputPath, 'cm-motion3', folderstr, sprintf('motion3-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        delete(C);
        
        % Top view.
        view(2);
        C = surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
        export_fig(fullfile(outputPath, 'cm-motion2', folderstr, sprintf('motion2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        
        % Plot divergence.
        figure(1);
        cla;
        colormap default;
        hold on;
        trisurf(F, S{t}(:, 1), S{t}(:, 2), S{t}(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sdiv{t, p});
        view(2);
        adjust3dplot;
        colorbar;
        if(sdivmax{p} - sdivmin{p} > eps)
            caxis([sdivmin{p}, sdivmax{p}]);
        end
        export_fig(fullfile(outputPath, 'cm-div2', folderstr, sprintf('div2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        close all;
    end
end
close all;