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
function tests = surfgradnormsquaredTest
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
theta = pi/1e-4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create Fourier coefficients for spherical surface.
c = 3.5449077018;

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute metric tensor.
f = surfgradnormsquared(N, c, [el, az]);
verifyEqual(testCase, size(f), [n, 1]);

% Check if spherical metric.
verifyEqual(testCase, f, zeros(n, 1));

end

function result2Test(testCase)

% Create triangulation of unit sphere.
[~, V] = sphTriang(4);
n = size(V, 1);

% Set parameters.
N = 0:1;

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create Fourier coefficients for spherical surface.
c = [3.5449077018, 1, 0, 0]';

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute metric tensor.
f = surfgradnormsquared(N, c, [el, az]);
verifyEqual(testCase, size(f), [n, 1]);

v = sum(surfgrad(N, c, [el, az]).^2, 2);
verifyEqual(testCase, v, f, 'absTol', 1e-6);

end

function visualizeTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Set parameters.
N = 0:2;

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e-4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create Fourier coefficients for spherical surface.
c = [3.5449077018, 1, 0, 0, 0, 0, 0, 0, 0]';

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute surface synthesis.
[S, rho] = surfsynth(N, V, c);

% Compute metric tensor.
f = surfgradnormsquared(N, c, [el, az]);

figure;
hold on;
subplot(1, 3, 1);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), rho, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('rho');
daspect([1, 1, 1]);
view(3);
subplot(1, 3, 2);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), f, 'EdgeColor', 'none', 'FaceColor', 'interp');
title('squared norm of surf gradient of rho');
daspect([1, 1, 1]);
view(3);
subplot(1, 3, 3);
trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'EdgeColor', 'none', 'FaceColor', 'interp');
title('surface');
daspect([1, 1, 1]);
view(3);

end