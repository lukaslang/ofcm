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
function test_suite = surfchristoffelTest
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

% Set radius.
r = 100;

% Create Fourier coefficients for spherical surface.
c = r*3.5449077018;

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute surface synthesis.
[~, rho] = surfsynth(N, V, c);

% Compute Christoffel symbols.
[c111, c112, c121, c122, c221, c222] = surfchristoffel(N, c, rho, [el, az]);
assertEqual(size(c111), [n, 1]);
assertEqual(size(c112), [n, 1]);
assertEqual(size(c121), [n, 1]);
assertEqual(size(c122), [n, 1]);
assertEqual(size(c221), [n, 1]);
assertEqual(size(c222), [n, 1]);

assertAlmostEqual(c111, zeros(n, 1));
assertAlmostEqual(c121, zeros(n, 1));
assertAlmostEqual(c221, -sin(el) .* cos(el));

assertAlmostEqual(c112, zeros(n, 1));
assertAlmostEqual(c122, cot(el));
assertAlmostEqual(c222, zeros(n, 1));

end

function visualizeTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Set parameters.
N = 0:1;

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create Fourier coefficients for spherical surface.
c = [3.5449077018, 0, 0, 0]';

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute surface synthesis.
[S, rho] = surfsynth(N, V, c);

% Compute Christoffel symbols.
[c111, c112, c121, c122, c221, c222] = surfchristoffel(N, c, rho, [el, az]);

figure;
hold on;
subplot(1, 4, 1);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  c111, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('c111');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 2);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  c121, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('c121');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 3);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  c121, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('c211');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 4);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  c221, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('c221');
daspect([1, 1, 1]);
view(3);

figure;
hold on;
subplot(1, 4, 1);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  c112, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('c112');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 2);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  c122, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('c122');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 3);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  c122, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('c212');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 4);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  c222, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('c222');
daspect([1, 1, 1]);
view(3);

end