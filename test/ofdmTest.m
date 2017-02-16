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
function tests = ofdmTest
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
k = 3;
h = 0.95;

% Define degree of integration.
deg = 100;

% Create sphere-like surface.
Ns = 0;
cs = 3.5449077018;

% Create integration points and quadrature rule for spherical cap.
[xi, w] = gausslegendre(deg, pi/2);
n = size(xi, 1);

% Triangulate of the upper unit hemi-sphere for placement of basis functions.
[~, X] = halfsphTriang(ref);

% Create data.
gradfx = zeros(n, 3);
dtfx = zeros(n, 1);

% Set regularisation parameter.
alpha = 1;

% Compute coefficients for optical flow.
ofc = ofdm(Ns, cs, X, k, h, xi, w, gradfx, dtfx, alpha);
verifyEqual(testCase, size(ofc), [2*size(X, 1), 1]);
verifyEqual(testCase, ofc, zeros(2*size(X, 1), 1));

end

function visualizeTest(testCase)

% Set subdivision parameter (number of basis functions is approx. 10*4^n).
ref = 4;

% Set parameters for basis function.
k = 4;
h = 0.95;

% Define degree of integration.
deg = 100;

% Create sphere.
Ns = 0;
cs = 3.5449077018;

% Create integration points and quadrature rule for spherical cap.
[xi, w] = gausslegendre(deg, pi/2);

% Compute evaluation points.
Y = sphcoord(xi, eye(3));

% Create rotation matrix.
theta = pi/9;
T = [   cos(theta),     sin(theta), 0;
       -sin(theta),     cos(theta), 0;
                 0,              0, 1];

% Create two images.
Ynj = spharm(1, Y);
f1 = Ynj(:, 1);
Ynj = spharm(1, Y*T);
f2 = Ynj(:, 1);
dtfx = f2 - f1;

N = 0:1;
c = [0, 1, 0, 0]';
gradfx = surfgrad(N, c, xi);

% Triangulate of the upper unit hemi-sphere for placement of basis functions.
[~, X] = halfsphTriang(ref);

% Set regularisation parameter.
alpha = 1;

% Compute coefficients for optical flow.
ofc = ofdm(Ns, cs, X, k, h, xi, w, gradfx, dtfx, alpha);

% Create triangulation for visualization purpose.
[F, V] = halfsphTriang(5);
[S, ~] = surfsynth(Ns, V, cs);

% Evaluate data at vertices.
Ynj = spharm(1, V);
fd{1} = Ynj(:, 1);
Ynj = spharm(1, V*T);
fd{2} = Ynj(:, 1);

% Find midpoints of faces on sphere.
TR = TriRep(F, V);
IC = normalise(TR.incenters);
[ICS, ~] = surfsynth(Ns, IC, cs);

% Evaluate basis functions at vertices.
[bfc1, bfc2] = vbasiscomp(k, h, X, IC);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(IC(:, 1), IC(:, 2), IC(:, 3));
el = pi/2 - el;

% Compute tangent basis.
[d1, d2] = surftanbasis(Ns, cs, [el, az]);

% Compute pushforward of basis functions.
v = bsxfun(@times, (bfc1')*ofc, d1) + bsxfun(@times, (bfc2')*ofc, d2);
clear bfc1;
clear bfc2;

% Visualize surface.
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
quiver3(ICS(:, 1), ICS(:, 2), ICS(:, 3), v(:, 1), v(:, 2), v(:, 3), 1, 'r');

figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{2}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);

% Project and scale flow.
up = projecttoplane(v);

% Compute colour space scaling.
nmax = max(sqrt(sum(up.^2, 2)));

% Compute colour of projection.
col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

% Visualize flow.
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', col, 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

% Plot colourwheel.
figure;
cw = colourWheel;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(2);

end

function visualize2Test(testCase)

% Set subdivision parameter (number of basis functions is approx. 10*4^n).
ref = 3;

% Set parameters for basis function.
k = 4;
h = 0.9;

% Define degree of integration.
deg = 300;

% Create sphere.
Ns = 0;
cs = 3.5449077018;

% Create integration points and quadrature rule for spherical cap.
[xi, w] = gausslegendre(deg, pi/2);

% Compute evaluation points.
Y = sphcoord(xi, eye(3));

% Create two images.
kd = 3;
hd = 0.99;
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

% Set regularisation parameter.
alpha = 1;

% Compute coefficients for optical flow.
ofc = ofdm(Ns, cs, X, k, h, xi, w, gradfx, dtfx, alpha);

% Create triangulation for visualization purpose.
[F, V] = halfsphTriang(5);
[S, ~] = surfsynth(Ns, V, cs);

% Evaluate data at vertices.
fd{1} = basisfun(kd, hd, sphcoord([0, 0], eye(3)), V);
fd{2} = basisfun(kd, hd, sphcoord([pi/1e2, 0], eye(3)), V);

% Find midpoints of faces on sphere.
TR = TriRep(F, V);
IC = normalise(TR.incenters);
[ICS, ~] = surfsynth(Ns, IC, cs);

% Evaluate basis functions at vertices.
[bfc1, bfc2] = vbasiscomp(k, h, X, IC);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(IC(:, 1), IC(:, 2), IC(:, 3));
el = pi/2 - el;

% Compute tangent basis.
[d1, d2] = surftanbasis(Ns, cs, [el, az]);

% Compute pushforward of basis functions.
v = bsxfun(@times, (bfc1')*ofc, d1) + bsxfun(@times, (bfc2')*ofc, d2);
clear bfc1;
clear bfc2;

% Visualize grad.
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
[Sy, ~] = surfsynth(Ns, Y, cs);
quiver3(Sy(:, 1), Sy(:, 2), Sy(:, 3), gradfx(:, 1), gradfx(:, 2), gradfx(:, 3), 1, 'r');

% Visualize flow.
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
quiver3(ICS(:, 1), ICS(:, 2), ICS(:, 3), v(:, 1), v(:, 2), v(:, 3), 0, 'r');

figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{2}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);

% Project and scale flow.
up = projecttoplane(v);

% Compute colour space scaling.
nmax = max(sqrt(sum(up.^2, 2)));

% Compute colour of projection.
col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

% Visualize flow.
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', col, 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

% Plot colourwheel.
figure;
cw = colourWheel;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(2);

end

function visualize3Test(testCase)

% Set subdivision parameter (number of basis functions is approx. 10*4^n).
ref = 3;

% Set parameters for basis function.
k = 4;
h = 0.9;

% Define degree of integration.
deg = 300;

% Create sphere.
Ns = 0;
cs = 3.5449077018;

% Create integration points and quadrature rule for spherical cap.
[xi, w] = gausslegendre(deg, pi/2);

% Compute evaluation points.
Y = sphcoord(xi, eye(3));

% Create two images.
kd = 3;
hd = 0.99;
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

% Set regularisation parameter.
alpha = 1;
beta = 1;

% Compute coefficients for optical flow.
% Compute optimality conditions.
[~, A, D, E, b] = optcond(Ns, cs, X, k, h, xi, w, gradfx, dtfx, fd{1});

% Solve linear system.
[ofc, ~] = solvesystem(A + alpha * D + beta * E, b, 1e-6, 1000);

% Create triangulation for visualization purpose.
[F, V] = halfsphTriang(5);
[S, ~] = surfsynth(Ns, V, cs);

% Evaluate data at vertices.
fd{1} = basisfun(kd, hd, sphcoord([0, 0], eye(3)), V);
fd{2} = basisfun(kd, hd, sphcoord([pi/1e2, 0], eye(3)), V);

% Find midpoints of faces on sphere.
TR = TriRep(F, V);
IC = normalise(TR.incenters);
[ICS, ~] = surfsynth(Ns, IC, cs);

% Evaluate basis functions at vertices.
[bfc1, bfc2] = vbasiscomp(k, h, X, IC);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(IC(:, 1), IC(:, 2), IC(:, 3));
el = pi/2 - el;

% Compute tangent basis.
[d1, d2] = surftanbasis(Ns, cs, [el, az]);

% Compute pushforward of basis functions.
v = bsxfun(@times, (bfc1')*ofc, d1) + bsxfun(@times, (bfc2')*ofc, d2);
clear bfc1;
clear bfc2;

% Visualize grad.
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
[Sy, ~] = surfsynth(Ns, Y, cs);
quiver3(Sy(:, 1), Sy(:, 2), Sy(:, 3), gradfx(:, 1), gradfx(:, 2), gradfx(:, 3), 1, 'r');

% Visualize flow.
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
quiver3(ICS(:, 1), ICS(:, 2), ICS(:, 3), v(:, 1), v(:, 2), v(:, 3), 0, 'r');

figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{1}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), fd{2}, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
view(3);

% Project and scale flow.
up = projecttoplane(v);

% Compute colour space scaling.
nmax = max(sqrt(sum(up.^2, 2)));

% Compute colour of projection.
col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

% Visualize flow.
figure;
hold on;
trisurf(F, S(:, 1), S(:, 2), S(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', col, 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

% Plot colourwheel.
figure;
cw = colourWheel;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(2);

end