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
function tests = cmTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Set subdivision parameter (number of basis functions is approx. 10*4^n).
ref = 3;

% Set parameters for basis function.
k = 4;
h = 0.9;

% Define degree of integration.
deg = 150;

% Create sphere.
Ns = 0;
cs = 3.5449077018;

% Create integration points and quadrature rule for spherical cap.
[xi, w] = gausslegendre(deg, pi/2);

% Compute evaluation points.
Y = sphcoord(xi, eye(3));

% Create two images.
kd = 3;
hd = 0.95;
fd{1} = basisfun(kd, hd, sphcoord([0, 0], eye(3)), Y)';
fd{2} = basisfun(kd, hd, sphcoord([pi/1e2, 0], eye(3)), Y)';
dtfx = fd{2} - fd{1};

% Compute gradient of first image.
%u = vbasisfun(kd, hd, sphcoord([0, 0], eye(3)), Y);
% Compute components of basis functions.
[dy1, dy2] = vbasiscomp(kd, hd, sphcoord([0, 0], eye(3)), Y);
% Compute coordinate basis at evaluation points.
[d1, d2] = surftanbasis(Ns, cs, xi);
% Compute gradient of basis function and orthogonal to gradient.
u = bsxfun(@times, full(dy1), permute(d1, [3, 1, 2])) + bsxfun(@times, full(dy2), permute(d2, [3, 1, 2]));
gradfx = squeeze(u(1, :, :));

% Triangulate of the upper unit hemi-sphere for placement of basis functions.
[~, X] = halfsphTriang(ref);

% Compute optimality conditions.
[~, A, D, E, G, b] = optcondcm(Ns, cs, cs, X, k, h, xi, w, gradfx, dtfx, fd{1}, fd{1}, 1024^3);
[~, ~, Acm, Dcm, Ecm, Gcm, ~, bcm] = optcondofcm(Ns, cs, cs, X, k, h, xi, w, gradfx, dtfx, fd{1}, fd{1}, 1024^3);
verifyEqual(testCase, Acm, A, 'absTol', 1e-15);
verifyEqual(testCase, Dcm, D, 'absTol', 1e-15);
verifyEqual(testCase, Ecm, E, 'absTol', 1e-15);
verifyEqual(testCase, Gcm, G, 'absTol', 1e-15);
verifyEqual(testCase, bcm, b, 'absTol', 1e-15);

% Compute optimality conditions.
[~, A, D, E, b] = optcond(Ns, cs, X, k, h, xi, w, gradfx, dtfx, fd{1}, 1024^3);
[~, Aof, ~, Dof, Eof, ~, bof, ~] = optcondofcm(Ns, cs, cs, X, k, h, xi, w, gradfx, dtfx, fd{1}, fd{1}, 1024^3);
verifyEqual(testCase, Aof, A, 'absTol', 1e-15);
verifyEqual(testCase, Dof, D, 'absTol', 1e-15);
verifyEqual(testCase, Eof, E, 'absTol', 1e-15);
verifyEqual(testCase, bof, b, 'absTol', 1e-15);

end