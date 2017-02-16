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
function tests = spharmderivnTest
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
Ns = 0:5;

% Create evaluation points.
[~, y] = sphTriang(4);
n = size(y, 1);

% Compute coordinates.
[az, el, ~] = cart2sph(y(:, 1), y(:, 2), y(:, 3));
el = pi/2 - el;

% Create basis functions.
[d1, d2] = spharmderivn(Ns, [el, az]);

m = Ns(end)^2 + 2*Ns(end) - Ns(1)^2 + 1;
verifyEqual(testCase, size(d1), [n, m]);
verifyEqual(testCase, size(d2), [n, m]);

end