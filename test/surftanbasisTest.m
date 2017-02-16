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
function tests = surftanbasisTest
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

% Compute tangent basis.
[d1, d2] = surftanbasis(N, c, [el, az]);
verifyEqual(testCase, size(d1), [n, 3]);
verifyEqual(testCase, size(d2), [n, 3]);

% Check if orthogonal to points.
verifyEqual(testCase, dot(d1, V, 2), zeros(n, 1), 'absTol', 1e-15);
verifyEqual(testCase, dot(d2, V, 2), zeros(n, 1), 'absTol', 1e-15);

% Compare to tangent basis of sphere.
[d1s, d2s] = sphtanbasis([el, az], eye(3));
verifyEqual(testCase, d1, d1s, 'absTol', 1e-10);
verifyEqual(testCase, d2, d2s, 'absTol', 1e-10);

end