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
function test_suite = basisfunTest
    initTestSuite;
end

function resultTest

% Set parameters.
k = 3;
h = 0;

% Create point x.
x = [0, 0, 1];
m = size(x, 1);

% Create evaluation points.
[~, y] = sphTriang(4);
n = size(y, 1);

% Create basis functions.
b = basisfun(k, h, x, y);
assertEqual(size(b), [m, n]);

end

function visualizeTest

% Set parameters.
k = 3;
h = 0.75;

% Create point x.
x = [0, 0, 1];
m = size(x, 1);

% Create evaluation points.
[F, V] = sphTriang(4);
n = size(V, 1);

% Create basis functions.
b = basisfun(k, h, x, V);
assertEqual(size(b), [m, n]);

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), b, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);

end