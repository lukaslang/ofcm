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
function tests = aligncentersTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Create centers.
S = cell(2, 1);
[X, Y, Z] = sphere(5);
S{1} = [X(:), Y(:), Z(:)];
[X, Y, Z] = sphere(6);
S{2} = [X(:), Y(:), Z(:)];

% Align.
[C, sc, sr] = aligncenters(S, 1.5);
verifyEqual(testCase, C{1}, S{1}, 'absTol', 1e-4);
verifyEqual(testCase, C{2}, S{2}, 'absTol', 1e-4);
verifyEqual(testCase, sc, [0, 0, 0], 'absTol', 1e-4);
verifyEqual(testCase, sr, 1, 'absTol', 1e-4);

end