% Copyright 2015 Lukas Lang
%
% This file is part of OFDM.
%
%    OFDM is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFDM is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFDM.  If not, see <http://www.gnu.org/licenses/>.
function test_suite = surfinvmetricTest
    initTestSuite;
end

function resultTest

% Create triangulation of unit sphere.
[~, V] = sphTriang(4);
n = size(V, 1);

% Set parameters.
N = 0;

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create Fourier coefficients for spherical surface.
c = 3.5449077018;

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute surface synthesis.
[~, rho] = surfsynth(N, V, c);

% Compute detg.
[g11, g12, g21, g22] = surfmetric(N, c, rho, [el, az]);
detg = g11 .* g22 - g12 .* g21;

% Compute inverse metric tensor.
[ginv11, ginv12, ginv21, ginv22] = surfinvmetric(N, c, rho, [el, az], detg);
assertEqual(size(ginv11), [n, 1]);
assertEqual(size(ginv12), [n, 1]);
assertEqual(size(ginv21), [n, 1]);
assertEqual(size(ginv22), [n, 1]);

% Check if spherical metric.
assertAlmostEqual(ginv12, zeros(n, 1));
assertAlmostEqual(ginv21, zeros(n, 1));
assertAlmostEqual(ginv11, ones(n, 1), 1e-6);
assertAlmostEqual(ginv22, 1 ./ sin(el).^2, 1e-6);

end

function visualizeTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Set parameters.
N = 0;

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create Fourier coefficients for spherical surface.
c = 3.5449077018;

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute surface synthesis.
[S, rho] = surfsynth(N, V, c);

% Compute detg.
[g11, g12, g21, g22] = surfmetric(N, c, rho, [el, az]);
detg = g11 .* g22 - g12 .* g21;

% Compute inverse metric tensor.
[ginv11, ginv12, ginv21, ginv22] = surfinvmetric(N, c, rho, [el, az], detg);

figure;
hold on;
subplot(1, 4, 1);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  ginv11, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('ginv11');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 2);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  ginv12, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('ginv12');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 3);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  ginv21, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('ginv21');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 4);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  ginv22, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('ginv22');
daspect([1, 1, 1]);
view(3);

end