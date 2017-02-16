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
function tests = suppTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Create one vertex.
V = [1, 0, 0];
h = 0;

% Compute support matrix.
S = supp(V, h);
verifyEqual(testCase, size(S), [1, 1]);
verifyEqual(testCase, S, true);

end

function result2Test(testCase)

% Create one vertex.
V = [1, 0, 0;
     0, 1, 0];
h = 0.75;

% Compute support matrix.
S = supp(V, h);
verifyEqual(testCase, size(S), [2, 2]);
verifyEqual(testCase, S, logical(eye(2)));

end

function result3Test(testCase)

% Create one vertex.
V = [1, 0, 0;
     0, 1, 0];
h = cos(pi/4);

% Compute support matrix.
S = supp(V, h);
verifyEqual(testCase, size(S), [2, 2]);
verifyEqual(testCase, S, true(2));

end