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
function test_suite = sphtanbasisTest
    initTestSuite;
end

function resultTest

% Define points.
xi = [0, 0;
      pi/2, 0;
      pi, pi];

% Compute tangent basis.  
[d1, d2] = sphtanbasis(xi, eye(3));
assertEqual(size(d1), [size(xi, 1), 3]);
assertEqual(size(d2), [size(xi, 1), 3]);

end

function visualizeTest

% Create evaluation points.
[F, V] = sphTriang(5);

% Convert to coordinates.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;

% Compute points.  
[d1, d2] = sphtanbasis([el, az], eye(3));

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
quiver3(V(:, 1), V(:, 2), V(:, 3), d1(:, 1), d1(:, 2), d1(:, 3), 'r');
quiver3(V(:, 1), V(:, 2), V(:, 3), d2(:, 1), d2(:, 2), d2(:, 3), 'b');
view(3);

end

function visualize2Test

% Create evaluation points.
[F, V] = sphTriang(5);

% Define rotation around y-axis.
theta = pi/4;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];

% Find coordinates.
[az, el, ~] = cart2sph(V(:, 1), V(:, 2), V(:, 3));
el = pi/2 - el;
[el, az] = convertsphcoord([el, az], R');

% Compute points.  
[d1, d2] = sphtanbasis([el, az], R');

figure;
hold on;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);
quiver3(V(:, 1), V(:, 2), V(:, 3), d1(:, 1), d1(:, 2), d1(:, 3), 'r');
quiver3(V(:, 1), V(:, 2), V(:, 3), d2(:, 1), d2(:, 2), d2(:, 3), 'b');
view(3);

end