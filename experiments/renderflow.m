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
[F, V] = halfsphTriang(6);

% Find midpoints of faces on sphere.
TR = TriRep(F, V);
IC = normalise(TR.incenters);

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

% Initialise.
for p=1:length(c)
    sdivmin{p} = inf;
    sdivmax{p} = -inf;
    wmin{p} = inf;
    wmax{p} = -inf;
    wpmax{p} = -inf;
    Umax{p} = -inf;
    Upmax{p} = -inf;
end

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
        sdiv = bfcdiv'*c{p};
        
        % Compute flow.
        w = bsxfun(@times, full((bfc1')*c{p}), d1) + bsxfun(@times, full((bfc2')*c{p}), d2);
        
        % Compute norm of flow.
        wnorm = sqrt(sum(w.^2, 2));
        
        % Project and scale flow.
        wp = projecttoplane(w);

        % Compute colour space scaling.
        wpnorm = sqrt(sum(wp.^2, 2));
        
        % Recover total velocity.
        U = Vs + w;
        Unorm = sqrt(sum(U.^2, 2));
        
        % Project and scale flow.
        Up = projecttoplane(U);
        Upnorm = sqrt(sum(Up.^2, 2));
        
        % Update min/max values.
        sdivmin{p} = min(sdivmin, min(sdiv));
        sdivmax{p} = max(sdivmax, max(sdiv));
        wmin{p} = min(wmin, min(wnorm));
        wmax{p} = max(wmax, max(wnorm));
        wpmax{p} = max(wmax, max(wpnorm));
        Umax{p} = max(Umax, max(Unorm));
        Upmax{p} = max(Upmax, max(Upnorm));        
    end
end
clear sdiv;
clear w;
clear wp;
clear U;
clear Up;
clear wnorm;
clear wpnorm;
clear Unorm;
clear Upnorm;

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
    
    % Compute synthesis for vertices.
    [S, ~] = surfsynth(Ns, V, cs{t});

    % Evaluate data at vertices.
    fd = cell2mat(evaldata(f(t), scale, {S}, sc, bandwidth, layers));

    % Compute synthesis for midpoints.
    [ICS, ~] = surfsynth(Ns, IC, cs{t});
    
    % Compute tangent basis.
    [d1, d2] = surftanbasis(Ns, cs{t}, xi);
    
    % Compute surface divergence of basis functions.
    bfcdiv = surfdivmem(Ns, cs{t}, xi, k, h, X, mem);
    
    % Compute surface normals.
    N = surfnormals(Ns, cs{t}, xi);

    % Compute surface velocity.
    dtrho = surfsynth(Ns, IC, (cs{t+1} - cs{t}) / dt);
    Vs = bsxfun(@times, dtrho, IC);
    
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
        
        % Compute flow.
        w = bsxfun(@times, full((bfc1')*c{p}), d1) + bsxfun(@times, full((bfc2')*c{p}), d2);
        
        % Project and scale flow.
        wp = projecttoplane(w);

        % Recover total velocity.
        U = Vs + w;
        
        % Project and scale flow.
        Up = projecttoplane(U);
        
        % Optical flow vectors.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd, 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        adjust3dplot;
        quiver3(ICS(:, 1), ICS(:, 2), ICS(:, 3), w(:, 1), w(:, 2), w(:, 3), 0, 'r');
        export_fig(fullfile(outputPath, 'of-vec2', folderstr, sprintf('vec2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(wp(:, 1)/wpmax{p}, wp(:, 2)/wpmax{p}))) ./ 255;

        % Optical flow.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
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
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sqrt(sum(w.^2, 2)));
        view(2);
        adjust3dplot;
        colorbar;
        if(wmax{p} - wmin{p} > eps)
            caxis([wmin{p}, wmax{p}]);
        end
        export_fig(fullfile(outputPath, 'of-flowmagn2', folderstr, sprintf('flowmagn2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        close all;

        % Compute colour of projection.
        col = double(squeeze(computeColour(Up(:, 1)/Upmax{p}, Up(:, 2)/Upmax{p}))) ./ 255;
        
        % Total velocity.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
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
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sdiv);
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
clear w;
clear Up;

% Initialise.
for p=1:length(c)
    sdivmin{p} = inf;
    sdivmax{p} = -inf;
    umin{p} = inf;
    umax{p} = -inf;
    upmax{p} = -inf;
    Umax{p} = -inf;
    Upmax{p} = -inf;
end

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
        sdiv = bfcdiv'*c{p};
        
        % Compute flow.
        u = bsxfun(@times, full((bfc1')*c{p}), d1) + bsxfun(@times, full((bfc2')*c{p}), d2);
        
        % Compute norm of flow.
        unorm = sqrt(sum(u.^2, 2));
        
        % Project and scale flow.
        up = projecttoplane(u);

        % Compute colour space scaling.
        upnorm = sqrt(sum(up.^2, 2));
        
        % Recover total velocity.
        U = bsxfun(@times, Vsn, N) + u;
        Unorm = sqrt(sum(U.^2, 2));
        
        % Project and scale flow.
        Up = projecttoplane(U);
        Upnorm = sqrt(sum(Up.^2, 2));
        
        % Update min/max values.
        sdivmin{p} = min(sdivmin, min(sdiv));
        sdivmax{p} = max(sdivmax, max(sdiv));
        umin{p} = min(umin, min(unorm));
        umax{p} = max(umax, max(unorm));
        upmax{p} = max(upmax, max(upnorm));
        Umax{p} = max(Umax, max(Unorm));
        Upmax{p} = max(Upmax, max(Upnorm));       
    end
end
clear sdiv;
clear u;
clear up;
clear U;
clear Up;
clear unorm;
clear upnorm;
clear Unorm;
clear Upnorm;

% Open a file to save flow scaling.
fid = fopen(fullfile(outputPath, 'velocities-cm.txt'), 'w');

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

    % Compute synthesis for vertices.
    [S, ~] = surfsynth(Ns, V, cs{t});

    % Evaluate data at vertices.
    fd = cell2mat(evaldata(f(t), scale, {S}, sc, bandwidth, layers));

    % Compute synthesis for midpoints.
    [ICS, ~] = surfsynth(Ns, IC, cs{t});
    
    % Compute tangent basis.
    [d1, d2] = surftanbasis(Ns, cs{t}, xi);
    
    % Compute surface divergence of basis functions.
    bfcdiv = surfdivmem(Ns, cs{t}, xi, k, h, X, mem);
    
    % Compute surface normals.
    N = surfnormals(Ns, cs{t}, xi);

    % Compute surface velocity.
    dtrho = surfsynth(Ns, IC, (cs{t+1} - cs{t}) / dt);
    Vs = bsxfun(@times, dtrho, IC);
    
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
        
        % Compute divergence.
        sdiv = bfcdiv'*c{p};
        
        % Compute flow.
        u = bsxfun(@times, full((bfc1')*c{p}), d1) + bsxfun(@times, full((bfc2')*c{p}), d2);
        
        % Project and scale flow.
        up = projecttoplane(u);

        % Recover total velocity.
        U = bsxfun(@times, Vsn, N) + u;
        
        % Project and scale flow.
        Up = projecttoplane(U);
        
        % Optical flow vectors.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd, 'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        adjust3dplot;
        quiver3(ICS(:, 1), ICS(:, 2), ICS(:, 3), u(:, 1), u(:, 2), u(:, 3), 0, 'r');
        export_fig(fullfile(outputPath, 'cm-vec2', folderstr, sprintf('vec2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(up(:, 1)/upmax{p}, up(:, 2)/upmax{p}))) ./ 255;

        % Optical flow.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
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
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sqrt(sum(u.^2, 2)));
        view(2);
        adjust3dplot;
        colorbar;
        if(umax{p} - umin{p} > eps)
            caxis([umin{p}, umax{p}]);
        end
        export_fig(fullfile(outputPath, 'cm-flowmagn2', folderstr, sprintf('flowmagn2-%s-%s-%.3i.png', name, folderstr, frames(t))), '-png', quality, '-transparent', '-a1');
        close all;
        
        % Compute colour of projection.
        col = double(squeeze(computeColour(Up(:, 1)/Upmax{p}, Up(:, 2)/Upmax{p}))) ./ 255;

        % Total velocity.
        figure(1);
        cla;
        colormap(cmap);
        hold on;
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
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
        trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', sdiv);
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