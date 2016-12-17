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
function test_suite = polarcoordTest
    initTestSuite;
end

function resultTest

% Create ONB.
e = eye(3);

% Define points.
xi = [1, 0;
      0, 0;
     -1, pi];

% Compute points.  
x = polarcoord(xi, e);
assertEqual(size(x), [size(xi, 1), 3]);
assertAlmostEqual(x, [0, 0, 1; 1, 0, 0; 0, 0, -1]);

% Check norm equals to one.
assertAlmostEqual(sum(x.^2, 2), ones(3, 1));

end

function result2Test

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
assertEqual(size(x), [size(xi, 1), 3]);
assertAlmostEqual(x, [-1, 0, 0; 0, 0, 1; 1, 0, 0]);

% Check norm equals to one.
assertAlmostEqual(sum(x.^2, 2), ones(3, 1));

end