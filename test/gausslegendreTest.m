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
function tests = gausslegendreTest
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
deg = 100;

% Compute number of points.
N = (floor(deg/2) + 1)*(deg + 1);

[xi, w] = gausslegendre(deg, pi/2);
verifyEqual(testCase, size(xi), [N, 2]);
verifyTrue(testCase, isvector(w));
verifyEqual(testCase, length(xi), N);
verifyEqual(testCase, length(w), N);

end

function integrateConstantFunctionTest(testCase)

% Set parameters.
deg = 100;

[~, w] = gausslegendre(deg, pi/2);

% Integrate constant function.
verifyEqual(testCase, sum(w), 2*pi, 'absTol', 1e-6);

end

function integrateConstantSphericalHarmonicTest(testCase)

% Set parameters.
deg = 100;
N = 0;

% Create quadrature for half-sphere.
[xi, w] = gausslegendre(deg, pi/2);

% Generate scalar spherical harmonics.
Y = spharmp(N, xi(:, 2), cos(xi(:, 1)));

% Integrate constant function.
verifyEqual(testCase, sum(Y.^2'*w), 0.5, 'absTol', 1e-6);

end

function integrateSphericalHarmonicTest(testCase)

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
    verifyEqual(testCase, v, 1, 'absTol', 1e-6);
end

end