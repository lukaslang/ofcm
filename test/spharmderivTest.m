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
function tests = spharmderivTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Set parameters.
Ns = 3;

% Create evaluation points.
[~, y] = sphTriang(4);
n = size(y, 1);

% Compute coordinates.
[az, el, ~] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
el = pi/2 - el;

% Create basis functions.
[d1, d2] = spharmderiv(Ns, [el, az]);
verifyEqual(testCase, size(d1), [n, 2*Ns + 1]);
verifyEqual(testCase, size(d2), [n, 2*Ns + 1]);

end

function visualizeTest(testCase)

% Set parameters.
N = 1;

% Create evaluation points.
[F, V] = sphTriang(4);

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/5e2;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Create basis functions.
Ynj = spharm(N, V);

% Compute coordinates.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

[d1, d2] = spharmderiv(N, [el, az]);

% Pick order.
k = 2;

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), Ynj(:, k), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);

figure;
hold on;
subplot(1, 2, 1);
trisurf(F, V(:, 1), V(:, 2), V(:, 3),  d1(:, k), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
subplot(1, 2, 2);
trisurf(F, V(:, 1), V(:, 2), V(:, 3),  d2(:, k), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);

end