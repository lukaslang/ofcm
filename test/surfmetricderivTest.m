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
function tests = surfmetricderivTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

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

% Compute derivatives of metric tensor.
[d1g11, d1g12, d1g21, d1g22, d2g11, d2g12, d2g21, d2g22] = surfmetricderiv(N, c, rho, [el, az]);
verifyEqual(testCase, size(d1g11), [n, 1]);
verifyEqual(testCase, size(d1g12), [n, 1]);
verifyEqual(testCase, size(d1g21), [n, 1]);
verifyEqual(testCase, size(d1g22), [n, 1]);
verifyEqual(testCase, size(d2g11), [n, 1]);
verifyEqual(testCase, size(d2g12), [n, 1]);
verifyEqual(testCase, size(d2g21), [n, 1]);
verifyEqual(testCase, size(d2g22), [n, 1]);

verifyEqual(testCase, d1g12, zeros(n, 1));
verifyEqual(testCase, d1g21, zeros(n, 1));
verifyEqual(testCase, d1g11, zeros(n, 1));
verifyEqual(testCase, d1g22, sin(2*el), 'absTol', 1e-10);
verifyEqual(testCase, d2g12, zeros(n, 1));
verifyEqual(testCase, d2g21, zeros(n, 1));
verifyEqual(testCase, d2g11, zeros(n, 1));
verifyEqual(testCase, d2g22, zeros(n, 1));

end

function visualizeTest(testCase)

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
c = [3.5449077018, 0, 0, 1]';

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute surface synthesis.
[S, rho] = surfsynth(N, V, c);

% Compute derivatives of metric tensor.
[d1g11, d1g12, d1g21, d1g22, d2g11, d2g12, d2g21, d2g22] = surfmetricderiv(N, c, rho, [el, az]);

figure;
hold on;
subplot(1, 4, 1);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  d1g11, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('d1g11');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 2);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  d1g12, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('d1g12');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 3);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  d1g21, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('d1g21');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 4);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  d1g22, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('d1g22');
daspect([1, 1, 1]);
view(3);

figure;
hold on;
subplot(1, 4, 1);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  d2g11, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('d2g11');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 2);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  d2g12, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('d2g12');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 3);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  d2g21, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('d2g21');
daspect([1, 1, 1]);
view(3);
subplot(1, 4, 4);
trisurf(F, S(:, 1), S(:, 2), S(:, 3),  d2g22, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('d2g22');
daspect([1, 1, 1]);
view(3);

end