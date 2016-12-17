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
function test_suite = convertsphcoordTest
    initTestSuite;
end

function resultTest

% Define points.
xi = [0, 0;
      pi/2, 0;
      pi, pi;
      pi/2, pi/2];

% Compute points.  
[el, az] = convertsphcoord(xi, eye(3));
assertAlmostEqual([el, az], xi);

end

function result2Test

% Define points.
xi = [0, 0;
      pi/2, 0;
      pi, pi;
      pi/2, pi/2];

% Create rotation matrix.
theta = pi/2;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];

% Compute points on sphere.
x = sphcoord(xi, eye(3));
assertAlmostEqual(x, [0, 0, 1; 1, 0, 0; 0, 0, -1; 0, 1, 0]);

% Compute points.
[el, az] = convertsphcoord(xi, R);
assertAlmostEqual([el, az], [pi/2, 0; pi, 0; pi/2, pi; pi/2, pi/2]);

% Check if points are equal.
y = sphcoord([el, az], R);
assertAlmostEqual(x, y);

end

function result3Test

% Create triangulation of unit sphere.
[~, V] = sphTriang(4);

% Compute coordinates of evaluation points.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;
xi = [el, az];

% Create rotation matrix.
theta = pi/4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];

% Compute points.
[el, az] = convertsphcoord(xi, R);

% Compute points on sphere.
y = sphcoord([el, az], R);

% Check if same as original.
assertAlmostEqual(V, y);

end