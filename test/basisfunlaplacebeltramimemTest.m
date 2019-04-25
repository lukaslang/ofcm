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
function tests = basisfunlaplacebeltramimemTest
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
k = 3;
h = 0.5;

% Create point x.
[~, x] = sphTriang(2);
m = size(x, 1);

% Create evaluation points.
[~, y] = sphTriang(4);
n = size(y, 1);

% Specify max. memory for matrix multiplication.
mem = 1*1024^3;

% Compute Laplace-Beltrami of basis functions.
b = basisfunlaplacebeltramimem(k, h, x, y, mem);
assertEqual(testCase, size(b), [m, n]);

end

function visualizeTest(testCase)

% Set parameters.
k = 3;
h = 0.75;

% Create point x.
x = [0, 0, 1];
m = size(x, 1);

% Create triangulation for visualization purpose.
[F, V] = sphTriang(4);
n = size(V, 1);

% Find points at north and south pole.
idx = V(:, 3) > 1-eps | V(:, 3) < -1+eps;

% Rotate away from singularities.
theta = pi/1e-4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
V(idx, :) = V(idx, :) * R;

% Specify max. memory for matrix multiplication.
mem = 1024^3;

% Compute basis functions.
b = full(basisfunmem(k, h, x, V, mem));
assertEqual(testCase, size(b), [m, n]);

% Compute Laplace-Beltrami of basis functions.
lb = full(basisfunlaplacebeltramimem(k, h, x, V, mem));
assertEqual(testCase, size(lb), [m, n]);

figure;
subplot(1, 2, 1);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), b, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
subplot(1, 2, 2);
trisurf(F, V(:, 1), V(:, 2), V(:, 3), lb, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);

end