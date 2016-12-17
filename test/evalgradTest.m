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
function test_suite = evalgradTest
    initTestSuite;
end

function resultTest

% Create evaluation points.
[~, V] = sphTriang(4);
n = size(V, 1);

% Create artificial data.
f = {ones(100, 100, 100); ones(100, 100, 100)};

% Set scale.
scale = [1, 1, 1] / 50;

% Set offset.
sc = [1, 1, 1];

% Evaluate.
gradfx = evalgrad(f, scale, {V; V}, {V; V}, sc, [0.8, 1], 10);
assertEqual(length(gradfx), 2);
assertEqual(size(gradfx{1}), [n, 3]);

end