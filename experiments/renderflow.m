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

% Define render quality (set to '-r600' for print quality).
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
[F, V] = halfsphTriang(7);
[S, ~] = cellfun(@(c) surfsynth(Ns, V, c), cs, 'UniformOutput', false);

% Evaluate data at vertices.
fd = evaldata(f, scale, S, sc, bandwidth, layers);

% Find midpoints of faces on sphere.
TR = TriRep(F, V);
IC = normalise(TR.incenters);
[ICS, ~] = cellfun(@(c) surfsynth(Ns, IC, c), cs, 'UniformOutput', false);

% Evaluate basis functions at vertices.
[bfc1, bfc2] = vbasiscompmem(k, h, X, IC, mem);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(IC(:, 1), IC(:, 2), IC(:, 3));
el = pi/2 - el;
xi = [el, az];

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

% Precompute.
for t=1:length(frames)-1
    % Load experiment.
    path = fullfile('results', name, timestamp1);
    file = fullfile(path, sprintf('%s-coeff-of-%.3i.mat', timestamp2, frames(t)));
    load(file);
    
    % Compute tangent basis.
    [d1, d2] = surftanbasis(Ns, cs{t}, xi);
    
    % Compute surface divergence of basis functions.
    bfcdiv = surfdivmem(Ns, cs{t}, xi, k, h, X, mem);
    
    % Compute surface normals.
    N = surfnormals(Ns, cs{t}, xi);

    % Compute surface velocity.
    dtrho = surfsynth(Ns, IC, (cs{t+1} - cs{t}) / dt);
    Vs = bsxfun(@times, dtrho, IC);

    for p=1:length(c)
        % Compute divergence.
        sdiv{t, p} = bfcdiv'*c{p};
        
        % Compute flow.
        w{t, p} = bsxfun(@times, full((bfc1')*c{p}), d1) + bsxfun(@times, full((bfc2')*c{p}), d2);
        
        % Compute norm of flow.
        wnorm{t, p} = sqrt(sum(w{t, p}.^2, 2));
        
        % Project and scale flow.
        wp{t, p} = projecttoplane(w{t, p});

        % Compute colour space scaling.
        wpnorm{t, p} = sqrt(sum(wp{t, p}.^2, 2));
        
        % Recover total velocity.
        U{t, p} = Vs + w{t, p};
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
    
    % Compute tangent basis.
    [d1, d2] = surftanbasis(Ns, cs{t}, xi);
    
    % Compute surface divergence of basis functions.
    bfcdiv = surfdivmem(Ns, cs{t}, xi, k, h, X, mem);
    
    % Compute surface normals.
    N = surfnormals(Ns, cs{t}, xi);

    % Compute surface velocity.
    dtrho = surfsynth(Ns, IC, (cs{t+1} - cs{t}) / dt);
    Vs = bsxfun(@times, dtrho, IC);
    
    % Compute scalar normal part of surface velocity.
    Vsn = dot(Vs, N);
    
    for p=1:length(c)
        % Compute divergence.
        sdiv{t, p} = bfcdiv'*c{p};
        
        % Compute flow.
        u{t, p} = bsxfun(@times, full((bfc1')*c{p}), d1) + bsxfun(@times, full((bfc2')*c{p}), d2);
        
        % Compute norm of flow.
        unorm{t, p} = sqrt(sum(u{t, p}.^2, 2));
        
        % Project and scale flow.
        up{t, p} = projecttoplane(u{t, p});

        % Compute colour space scaling.
        upnorm{t, p} = sqrt(sum(up{t, p}.^2, 2));
        
        % Recover total velocity.
        U{t, p} = bsxfun(@times, Vsn, N) + u{t, p};
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