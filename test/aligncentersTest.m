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
function test_suite = aligncentersTest
    initTestSuite;
end

function resultTest

% Create centers.
S = cell(2, 1);
[X, Y, Z] = sphere(5);
S{1} = [X(:), Y(:), Z(:)];
[X, Y, Z] = sphere(6);
S{2} = [X(:), Y(:), Z(:)];

% Align.
[C, sc, sr] = aligncenters(S, 1.5);
assertAlmostEqual(C{1}, S{1}, 1e-4);
assertAlmostEqual(C{2}, S{2}, 1e-4);
assertAlmostEqual(sc, [0, 0, 0], 1e-4);
assertAlmostEqual(sr, 1, 1e-4);

end