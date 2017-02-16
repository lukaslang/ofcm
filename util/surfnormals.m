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
function N = surfnormals(Ns, cs, xi)
%SURFNORMALS Computes surface normals.
%   
%   N = SURFNORMALS(Ns, cs, xi) takes degrees Ns and corresponding 
%   Fourier coefficients cs and evaluation points xi and returns the
%   surface normals.
%
%   Ns is a list os non-negative consecutive integers.
%   cs is a vector of Fourier coefficients.
%   xi is an n-by-2 matrix of coordinates [el, az], where el in [0, pi] and
%   az in [0, 2pi).
%
%   N is a n-by-3 matrices.

[d1, d2] = surftanbasis(Ns, cs, xi);
N = cross(d1, d2, 2);
N = bsxfun(@rdivide, N, sqrt(sum(N .^ 2, 2)));

end