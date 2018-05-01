function plotflow(F, S, ICS, fd1, fd2, v, cmap, plotflowstreamlines)
%PLOTFLOW Creates figures and plots results.

% Set face alpha for detail.
faceAlpha = 0.5;

% Set colormap for streamlines.
cmaps = 'summer';

% Create colourwheel.
cw = colourwheelbg;

% Project and scale flow.
vp = projecttoplane(v);

% Compute colour space scaling.
R = max(sqrt(sum(vp.^2, 2)));

% Compute colour of projection.
col = double(squeeze(computeColour(vp(:, 1)./R, vp(:, 2)./R))) ./ 255;

% Compute length of longest projected tangent vector.
lmax = max(sqrt(sum(v(:, 1:2).^2, 2)));

% Set parameters for streamline computation.
stepsize = 1/lmax;
maxit = 50;

% Set line with of arrows and streamlines.
lineWidth = 1;

% Plot colour-coded flow.
figure(1);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
text(400, -350, 50, sprintf('R = %.2f', R), 'Interpreter', 'tex', 'FontName', 'Helvetica', 'FontSize', 14, 'Rotation', 90, 'HorizontalAlignment', 'left');
view(3);
adjust3dplot;

% Top view.
figure(2);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', col);
surf(250:399, -399:-250, zeros(150, 150), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
text(250, -370, 0, sprintf('R = %.2f', R), 'Interpreter', 'tex', 'FontName', 'Helvetica', 'FontSize', 14, 'HorizontalAlignment', 'right');
adjust3dplot;
view(2);

% First cell division.
figure(3);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd1, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', faceAlpha);
quiver3(ICS(:, 1), ICS(:, 2), ICS(:, 3), v(:, 1), v(:, 2), v(:, 3), 0, 'k', 'LineWidth', lineWidth);
adjust3dplot;
axis([-130, -60, 115, 185]);
%set(gca, 'XTick', -200:20:200);
%set(gca, 'YTick', -200:20:200);
view(2);
axis off;

figure(4);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd2, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', faceAlpha);
adjust3dplot;
axis([-130, -60, 115, 185]);
%set(gca, 'XTick', -200:20:200);
%set(gca, 'YTick', -200:20:200);
view(2);
axis off;

% Set seed points for streamlines.
[xs, ys] = meshgrid(-130:4:-60, 115:4:185);

% Create fake grid to prevent cropping.
[X, Y] = meshgrid(-130:-60, 115:185);

% Plot streamlines for flow.
figure(7);
colormap(cmaps);
hold on;
surf(X, Y, zeros(size(X)), 'EdgeColor', 'none', 'FaceColor', 'w');
view(2);
adjust3dplot;
axis([-130, -60, 115, 185]);
%set(gca, 'XTick', -200:20:200);
%set(gca, 'YTick', -200:20:200);
streamlines2(ICS(:, 1:2), v(:, 1:2), [xs(:), ys(:)], stepsize, maxit, cmaps, lineWidth);
axis off;

% Second cell division.
figure(5);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd1, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', faceAlpha);
quiver3(ICS(:, 1), ICS(:, 2), ICS(:, 3), v(:, 1), v(:, 2), v(:, 3), 0, 'k', 'LineWidth', lineWidth);
adjust3dplot;
axis([115, 185, -135, -65]);
%set(gca, 'XTick', -200:20:200);
%set(gca, 'YTick', -200:20:200);
view(2);
axis off;

figure(6);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd2, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', faceAlpha);
adjust3dplot;
axis([115, 185, -135, -65]);
%set(gca, 'XTick', -200:20:200);
%set(gca, 'YTick', -200:20:200);
view(2);
axis off;

% Set seed points for streamlines.
[xs, ys] = meshgrid(115:4:185, -135:4:-65);

% Create fake grid to prevent cropping.
[X, Y] = meshgrid(115:185, -135:-65);

% Plot streamlines for flow.
figure(8);
colormap(cmaps);
hold on;
surf(X, Y, zeros(size(X)), 'EdgeColor', 'none', 'FaceColor', 'w');
view(2);
adjust3dplot;
axis([115, 185, -135, -65]);
%set(gca, 'XTick', -200:20:200);
%set(gca, 'YTick', -200:20:200);
streamlines2(ICS(:, 1:2), v(:, 1:2), [xs(:), ys(:)], stepsize, maxit, cmaps, lineWidth);
axis off;

if(plotflowstreamlines)
    % Set parameters for streamline computation.
    stepsize = 15/(4*lmax);
    maxit = 50;
    
    % Set seed points for streamlines.
    [xs, ys] = meshgrid(-400:15:400, -400:15:400);

    % Plot streamlines for flow.
    figure(9);
    colormap(cmaps);
    hold on;
    view(2);
    adjust3dplot;
    streamlines2(ICS(:, 1:2), v(:, 1:2), [xs(:), ys(:)], stepsize, maxit, cmaps, lineWidth);
end

end