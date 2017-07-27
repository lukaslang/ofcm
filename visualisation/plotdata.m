function plotdata(F, V, S, Xm, Ym, Zm, f, fd, rho, cmap, scale, sc, bandwidth)
%PLOTDATA Creates figures and plots results.

% Plot function rho on the unit sphere.
figure(1);
daspect([1, 1, 1]);
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), rho, 'EdgeColor', 'none');
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
caxis([300, 400]);
cbar = findobj(figure(1), 'tag', 'Colorbar');
set(cbar, 'YTick', 280:20:460);
set(cbar, 'TickLength', 0.02, 'YColor', [0, 0, 0]);
view(2);

% Raw data.
figure(2);
colormap(cmap);
hold on;
vol3d('cdata', f, 'XData', scale(1) * [0, size(f, 1)] - sc(1), 'YData', scale(2) * [0, size(f, 2)] - sc(2), 'ZData', scale(3) * [0, size(f, 3)] - sc(3));
view(3);
adjust3dplot;

% Raw data and surface with data.
figure(3);
colormap(cmap);
hold on;
vol3d('cdata', f, 'XData', scale(1) * [0, size(f, 1)] - sc(1), 'YData', scale(2) * [0, size(f, 2)] - sc(2), 'ZData', scale(3) * [0, size(f, 3)] - sc(3));
view(3);
adjust3dplot;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd, 'EdgeColor', 'none', 'FaceColor', 'interp');

% Add grid.
figure(4);
colormap(cmap);
hold on;
vol3d('cdata', f, 'XData', scale(1) * [0, size(f, 1)] - sc(1), 'YData', scale(2) * [0, size(f, 2)] - sc(2), 'ZData', scale(3) * [0, size(f, 3)] - sc(3));
view(3);
adjust3dplot;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd, 'EdgeColor', 'none', 'FaceColor', 'interp');
surf(Xm, Ym, Zm, ones(size(Xm)), 'EdgeColor', [0.7, 0.7, 0.7], 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeAlpha', 1);

% Cross-section of narrow band.
figure(5);
colormap(cmap);
hold on;
vol3d('cdata', f, 'XData', scale(1) * [0, size(f, 1)] - sc(1), 'YData', scale(2) * [0, size(f, 2)] - sc(2), 'ZData', scale(3) * [0, size(f, 3)] - sc(3));
view([-90, 0]);
adjust3dplot;
set(gca, 'XLim', [0, 50]);

% Define narrow band around surface.
Sinner = bandwidth(1) * S;
Souter = bandwidth(2) * S;
trisurf(F, Sinner(:, 1), Sinner(:, 2), Sinner(:, 3), zeros(size(fd, 1), 1), 'EdgeColor', 'green', 'LineWidth', 1.5);
trisurf(F, Souter(:, 1), Souter(:, 2), Souter(:, 3), zeros(size(fd, 1), 1), 'EdgeColor', 'red', 'LineWidth', 1.5);
trisurf(F, S(:, 1), S(:, 2), S(:, 3), zeros(size(fd, 1), 1), 'EdgeColor', [0.7, 0.7, 0.7], 'LineWidth', 1.5);

% Surface data.
figure(6);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd, 'EdgeColor', 'none', 'FaceColor', 'interp');
view(3);
adjust3dplot;

% Top view.
figure(7);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd, 'EdgeColor', 'none', 'FaceColor', 'interp');
adjust3dplot;
view(2);

% Create a plot with grid.
figure(8);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd, 'EdgeColor', 'none', 'FaceColor', 'interp');
adjust3dplot;
view(3);
surf(Xm, Ym, Zm, ones(size(Xm)), 'EdgeColor', [0.7, 0.7, 0.7], 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeAlpha', 1);

% Top view.
figure(9);
colormap(cmap);
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd, 'EdgeColor', 'none', 'FaceColor', 'interp');
adjust3dplot;
surf(Xm, Ym, Zm, ones(size(Xm)), 'EdgeColor', [0.7, 0.7, 0.7], 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeAlpha', 1);
view(2);

end