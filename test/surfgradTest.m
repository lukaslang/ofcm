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
function tests = surfgradTest
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

% Compute metric tensor.
f = surfgrad(N, c, [el, az]);
verifyEqual(testCase, size(f), [n, 3]);

% Check if zero.
verifyEqual(testCase, f, zeros(n, 3));

end

function visualizeTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
n = size(V, 1);

% Set parameters.
N = 0:3;

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create Fourier coefficients for spherical surface.
c = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]';

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute surface synthesis.
[~, rho] = surfsynth(N, V, c);

% Compute surface gradient.
f = surfgrad(N, c, [el, az]);

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), rho, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
quiver3(V(:, 1), V(:, 2), V(:, 3), f(:, 1), f(:, 2), f(:, 3), 'r');

% Check if tangential.
verifyEqual(testCase, dot(f, V, 2), zeros(n, 1), 'absTol', 1e-15);

end