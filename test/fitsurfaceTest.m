% Copyright 2013 Clemens Kirisits and Lukas Lang
%
% This file is part of OFD.
%
%    OFD is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFD is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFD.  If not, see <http://www.gnu.org/licenses/>.
function tests = fitsurfaceTest
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
N = 0:5;
s = 1;
beta0 = 1;
beta1 = 1;

% Fit surface to sphere.
[c, Y] = fitsurface(N, {V}, beta0, beta1, s);
verifyFalse(testCase, isempty(c));
verifyFalse(testCase, isempty(Y));

% Recover function at specified vertices.
rho = Y{1}*c{1};

% Compute coordinates of fitted surface at data points.
F = bsxfun(@times, V, rho);

% Check if points are still on unit sphere.
len = sqrt(sum(F .^ 2,2));
verifyEqual(testCase, len, ones(n, 1), 'absTol', 1e-12);

% Check if points equal.
verifyEqual(testCase, F, V, 'absTol', 1e-12);

end

function result2Test(testCase)

% Create triangulation of unit sphere.
[~, V] = sphTriang(4);

% Create data using one spherical harmonic of degree n and order m.
n = 5;
m = 3;
Y = spharm(n, V);
X = bsxfun(@times, V, 1 + Y(:, m));

% Set parameters.
N = 0:5;
s = 1;
beta0 = 1;
beta1 = 0;

% Fit surface.
[c, Y] = fitsurface(N, {X}, beta0, beta1, s);
verifyFalse(testCase, isempty(c));
verifyFalse(testCase, isempty(Y));

% Compute surface with old function.
[d, Z] = surffit(N, X, beta0, s);

% Check if points are almost equal.
verifyEqual(testCase, c{1}, d, 'absTol', 1e-6);
verifyEqual(testCase, Y{1}, Z, 'absTol', 1e-6);

end

function result3Test(testCase)

% Create triangulation of unit sphere.
[~, V] = sphTriang(4);

% Create data using one spherical harmonic of degree n and order m.
n = 5;
m = 3;
Y = spharm(n, V);
X = bsxfun(@times, V, 1 + Y(:, m));

% Set parameters.
N = 0:5;
s = 1;
beta0 = 1;
beta1 = 1;

% Fit surface.
[c, Y] = fitsurface(N, repmat({X}, 3, 1), beta0, beta1, s);
verifyFalse(testCase, isempty(c));
verifyFalse(testCase, isempty(Y));

% Compute surface with old function.
[d, Z] = surffit(N, X, beta0, s);

% Check if points are almost equal.
for k=1:3
    verifyEqual(testCase, c{k}, d, 'absTol', 1e-6);
    verifyEqual(testCase, Y{k}, Z, 'absTol', 1e-6);
end

end