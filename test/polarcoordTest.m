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
function tests = polarcoordTest
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
xi = [1, 0;
      0, 0;
     -1, pi];

% Compute points.  
x = polarcoord(xi, e);
verifyEqual(testCase, size(x), [size(xi, 1), 3]);
verifyEqual(testCase, x, [0, 0, 1; 1, 0, 0; 0, 0, -1]);

% Check norm equals to one.
verifyEqual(testCase, sum(x.^2, 2), ones(3, 1));

end

function result2Test(testCase)

% Create ONB.
e = [0, 0, 1;
     0, 1, 0;
    -1, 0, 0];

% Define points.
xi = [1, 0;
      0, 0;
     -1, pi];

% Compute points.  
x = polarcoord(xi, e);
verifyEqual(testCase, size(x), [size(xi, 1), 3]);
verifyEqual(testCase, x, [-1, 0, 0; 0, 0, 1; 1, 0, 0]);

% Check norm equals to one.
verifyEqual(testCase, sum(x.^2, 2), ones(3, 1));

end