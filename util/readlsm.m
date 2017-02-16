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
function [lsm, rng] = readlsm(file)
%READLSM Reads an LSM file.
%   
%   [lsm, rng] = readlsm(file) takes a file string and returns the actual
%   LSM data and a vector rng containing the frame numbers.
%
%   Note that this function uses tiffread.m version v3.31.
%   See http://www.cytosim.org/misc/index.html for further details.

lsm = tiffread(file);
rng = 1:lsm(1).lsm.DimensionTime;

end