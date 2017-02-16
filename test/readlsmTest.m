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
function tests = readlsmTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Specify dataset.
name = 'cxcr4aMO2_290112_Maximumintensityprojection';
dataFolder = fullfile(datapath, 'LSM 16.03.2012');
file = fullfile(dataFolder, strcat(name, '.lsm'));

% Read LSM file.
[lsm, rng] = readlsm(file);

verifyFalse(testCase, isempty(lsm));
verifyFalse(testCase, isempty(rng));
verifyEqual(testCase, rng, 1:151);

end