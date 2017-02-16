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
function [el, az] = convertsphcoord(xi, R)
%CONVERTSPHCOORD Converts spherical coordinates in the standard
%parametrisation to spherical coordinates in a rotated parametrisation.
%   
%   [el, az] = CONVERTSPHCOORD(xi, R) takes spherical coordinates xi and 
%   a rotation matrix R and returns coordinates for the same points on the 
%   2-sphere but in the coordinate system specified by R.
%
%   xi is a n-by-2 matrix with elements in [0, pi] and [0, 2*pi).
%   R is a 3-by-3 rotation matrix.
%   [el, az] is a n-by-2 matrix of coordinates.

% Compute points on sphere.
x = sphcoord(xi, eye(3));
    
% Rotate back.
y = x*R';

% Compute new coordinates.
el = acos(y(:, 3));
az = atan2(y(:, 2), y(:, 1));

end