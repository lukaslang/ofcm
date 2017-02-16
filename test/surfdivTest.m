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
function tests = surfdivTest
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
Ns = 0;

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e-4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Set radius.
r = 10;

% Create Fourier coefficients for spherical surface.
cs = r*3.5449077018;

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Set parameters.
k = 3;
h = 0;

% Create point x.
X = [0, 0, 1];

% Set coefficients
c = [0, 1]';

% Compute divergence.
v = surfdiv(Ns, cs, [el, az], k, h, X);

% Multiply with coefficients.
v = (v') * c;
verifyEqual(testCase, size(v), [n, 1]);
verifyEqual(testCase, v, zeros(n, 1));

end

function visualizeTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);

% Set parameters.
Ns = 0;

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e-4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Set radius.
r = 1;

% Create Fourier coefficients for spherical surface.
cs = r*3.5449077018;

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Set parameters.
k = 3;
h = 0.75;

% Create point x.
X = [0, 0, 1];

% Set coefficients
c = [1, 0]';

% Compute divergence.
v = surfdiv(Ns, cs, [el, az], k, h, X);

% Multiply with coefficients.
v = (v') * c;

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), v, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);

end