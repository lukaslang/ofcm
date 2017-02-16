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
function tests = sphcoordTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Create ONB.
e = eye(3);

% Define points.
xi = [0, 0;
      pi/2, 0;
      pi, pi];

% Compute points.  
x = sphcoord(xi, e);
verifyEqual(testCase, size(x), [size(xi, 1), 3]);
verifyEqual(testCase, x, [0, 0, 1; 1, 0, 0; 0, 0, -1], 'absTol', 1e-15);

% Check norm equals to one.
verifyEqual(testCase, sum(x.^2, 2), ones(3, 1), 'absTol', 1e-15);

end

function result2Test(testCase)

% Create ONB.
e = [0, 0, 1;
     0, 1, 0;
    -1, 0, 0];

% Define points.
xi = [0, 0;
      pi/2, 0;
      pi, pi];

% Compute points.  
x = sphcoord(xi, e);
verifyEqual(testCase, size(x), [size(xi, 1), 3]);
verifyEqual(testCase, x, [-1, 0, 0; 0, 0, 1; 1, 0, 0], 'absTol', 1e-15);

% Check norm equals to one.
verifyEqual(testCase, sum(x.^2, 2), ones(3, 1), 'absTol', 1e-15);

end

function result3Test(testCase)

% Create triangulation of unit sphere.
[~, V] = sphTriang(4);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Create rotation matrix.
theta = pi/4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];

% Compute points on sphere.
y = sphcoord([el, az], R);

% Check if same as rotated original.
verifyEqual(testCase, V*R, y, 'absTol', 1e-15);

end