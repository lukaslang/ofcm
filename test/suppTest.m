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
function test_suite = suppTest
    initTestSuite;
end

function resultTest

% Create one vertex.
V = [1, 0, 0];
h = 0;

% Compute support matrix.
S = supp(V, h);
assertEqual(size(S), [1, 1]);
assertEqual(S, true);

end

function result2Test

% Create one vertex.
V = [1, 0, 0;
     0, 1, 0];
h = 0.75;

% Compute support matrix.
S = supp(V, h);
assertEqual(size(S), [2, 2]);
assertEqual(S, logical(eye(2)));

end

function result3Test

% Create one vertex.
V = [1, 0, 0;
     0, 1, 0];
h = cos(pi/4);

% Compute support matrix.
S = supp(V, h);
assertEqual(size(S), [2, 2]);
assertEqual(S, true(2));

end