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
function [C, sc, sr] = aligncenters(C, sr)
%ALIGNCENTERS Returns approximate cell centres.
%   
%   [C, sc, r] = ALIGNCENTERS(C) takes a cell array of coordinates 
%   [X, Y, Z] of cell centers and aligns them around the origin. C again is
%   a cell array of same lenght as the input, sc the spherical center 
%   [x, y, z] of a fitted sphere, and sr its radius. sr is the initial
%   radius for the spherical fitting.

% Convert to matrix.
X = cell2mat(C);

% Initialize spherical fitting.
sc = mean(X);

% Fit sphere to all centers simultaneously.
[sc, sr] = spherefit(X, sc, sr);

% Subtract spherical center.
C = cellfun(@(x) bsxfun(@minus, x, sc), C, 'UniformOutput', false);

end