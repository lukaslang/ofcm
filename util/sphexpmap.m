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
function P = sphexpmap(V, v)
%SPHEXPMAP Computes the spherical exponential map.
%   
%   P = SPHEXPMAP(V, v, t) takes a matrix V of points on the sphere, a 
%   tangent vector field v at V, and computes the exponential map exp_{V}(v).
%
%   V is a matrix of size [n, 3].
%   v is a matrix of size [n, 3].
%
%   See https://math.stackexchange.com/questions/1923416/exponential-map-on-the-n-sphere

len = sqrt(sum(v.^2, 2));
idx = len > 1e-10;
P = V;
P(idx, :) = cos(len(idx)) .* V(idx, :) + sin(len(idx)) .* v(idx, :) ./ len(idx);

end