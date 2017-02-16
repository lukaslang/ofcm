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
function f = lsmframe(lsm, frame, band)
%LSMFRAME Returns a specified frame from LSM data.
%   
%   f = LSMFRAME(lsm, frame, band) takes lsm data and returns a frame from
%   a specified band.

% Get resolution.
x = lsm(1).lsm.DimensionX;
y = lsm(1).lsm.DimensionY;
z = lsm(1).lsm.DimensionZ;

% Compute offset.
offset = (frame - 1) * z;

% Check if multiple bands are present and get slices.
if(iscell(lsm(1).data))
    f = cell2mat(cellfun(@(a) a{1, band}, {lsm(offset + 1:offset + z).data}, 'UniformOutput', false));
else
    f = [lsm(offset + 1:offset + z).data];
end

% Reshape to frame size.
f = reshape(f, [x, y, z]);

end