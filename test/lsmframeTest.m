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
function test_suite = lsmframeTest
    initTestSuite;
end

function resultTest

% Create data.
d = randi(255, [512, 512, 10]);

% Create mock LSM data with a single frame.
for k=1:10
    lsm(k).lsm.DimensionX = 512;
    lsm(k).lsm.DimensionY = 512;
    lsm(k).lsm.DimensionZ = 10;
    lsm(k).data = d(:, :, k);
end

f = lsmframe(lsm, 1, 1);
assertEqual(f, d);

end

function multipleFramesTest

% Create data.
d = randi(255, [512, 512, 50]);

% Create mock LSM data with five frames.
for k=1:50
    lsm(k).lsm.DimensionX = 512;
    lsm(k).lsm.DimensionY = 512;
    lsm(k).lsm.DimensionZ = 10;
    lsm(k).data = d(:, :, k);
end

f = lsmframe(lsm, 2, 1);
assertEqual(f, d(:, :, 11:20));

end

function twoDimensionsTest

% Create data.
d = randi(255, [512, 512, 5]);

% Create mock LSM data with five frames.
for k=1:5
    lsm(k).lsm.DimensionX = 512;
    lsm(k).lsm.DimensionY = 512;
    lsm(k).lsm.DimensionZ = 1;
    lsm(k).data = d(:, :, k);
end

f = lsmframe(lsm, 3, 1);
assertEqual(f, d(:, :, 3));

end

function multipleFramesWithBandsTest

% Create data.
d1 = randi(255, [512, 512, 50]);
d2 = randi(255, [512, 512, 50]);

% Create mock LSM data with five frames and two bands.
for k=1:50
    lsm(k).lsm.DimensionX = 512;
    lsm(k).lsm.DimensionY = 512;
    lsm(k).lsm.DimensionZ = 10;
    lsm(k).data{1} = d1(:, :, k);
    lsm(k).data{2} = d2(:, :, k);
end

f = lsmframe(lsm, 2, 1);
assertEqual(f, d1(:, :, 11:20));

f = lsmframe(lsm, 3, 2);
assertEqual(f, d2(:, :, 21:30));

end