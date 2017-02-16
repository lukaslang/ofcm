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
function [f, scale] = loaddata(file, band, frames)
%LOADDATA Reads an LSM file and returns the data.
%   
%   f = LOADDATA(file, band, frames) takes a filename and reads lsm data 
%   and returns a cell array of volumetric images from a specified band.
%   frames is an array of frames to extract from the sequence.

% Load LSM file.
lsm = readlsm(file);

% Read LSM info.
lsminfo = lsm.lsm;

% Get scaling of data.
scale = [lsminfo.VoxelSizeX, lsminfo.VoxelSizeY, lsminfo.VoxelSizeZ];

% Create cell array.
nFrames = length(frames);
f = cell(nFrames, 1);

% Store each frame.
for k=1:nFrames
    f{k} = lsmframe(lsm, frames(k), band);
end

end