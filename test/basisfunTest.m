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
function tests = basisfunTest
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
h = 0;

% Create point x.
x = [0, 0, 1];
m = size(x, 1);

% Create evaluation points.
[~, y] = sphTriang(4);
n = size(y, 1);

% Create basis functions.
b = basisfun(k, h, x, y);
assertEqual(testCase, size(b), [m, n]);

end

function visualizeTest(testCase)

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
assertEqual(testCase, size(b), [m, n]);

figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), b, 'EdgeColor', 'none', 'FaceColor', 'interp');
daspect([1, 1, 1]);

end