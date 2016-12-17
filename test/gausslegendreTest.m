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
function test_suite = gausslegendreTest
    initTestSuite;
end

function resultTest

% Set parameters.
deg = 100;

% Compute number of points.
N = (floor(deg/2) + 1)*(deg + 1);

[xi, w] = gausslegendre(deg, pi/2);
assertEqual(size(xi), [N, 2]);
assertTrue(isvector(w));
assertEqual(length(xi), N);
assertEqual(length(w), N);

end

function integrateConstantFunctionTest

% Set parameters.
deg = 100;

[~, w] = gausslegendre(deg, pi/2);

% Integrate constant function.
assertAlmostEqual(sum(w), 2*pi);

end

function integrateConstantSphericalHarmonicTest

% Set parameters.
deg = 100;
N = 0;

% Create quadrature for half-sphere.
[xi, w] = gausslegendre(deg, pi/2);

% Generate scalar spherical harmonics.
Y = spharmp(N, xi(:, 2), cos(xi(:, 1)));

% Integrate constant function.
assertAlmostEqual(sum(Y.^2'*w), 0.5);

end

function integrateSphericalHarmonicTest

% Set parameters.
deg = 40;
N = 20;

% Create quadrature rule for almost whole sphere.
[xi, w] = gausslegendre(deg, pi-2*eps);

% Generate scalar spherical harmonics.
Y = spharmp(N, xi(:, 2), cos(xi(:, 1)));

for k=1:2*N+1
    % Integrate.
    v = (Y(:, k).^2)'*w;
    assertAlmostEqual(v, 1);
end

end