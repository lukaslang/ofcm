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
function tests = sphexpmapTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

V = [0, 0, 1];
v = 2*pi*[1, 0, 0];
P = sphexpmap(V, v);
verifyEqual(testCase, P, V, 'absTol', 1e-15);

V = [0, 0, 1];
v = 2*pi*[0, 1, 0];
P = sphexpmap(V, v);
verifyEqual(testCase, P, V, 'absTol', 1e-15);

V = [1, 0, 0];
v = pi*[1, 0, 0];
P = sphexpmap(V, v);
verifyEqual(testCase, P, -V, 'absTol', 1e-15);

V = sphcoord([pi/4, 0], eye(3));
v = pi*[0, 1, 0];
P = sphexpmap(V, v);
verifyEqual(testCase, P, -V, 'absTol', 1e-15);

end