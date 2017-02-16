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
function S = supp(V, h)
%SUPP Returns matrix of shared support.
%   
%   S = SUPP(V, h) takes a n-by-3 matrix V and a scaling parameter h and 
%   returns a symmetric n-by-n logical matrix S which indicates which
%   basis functions have shared support.

% Define dot product.
fun = @(x, y) sum(repmat(x, size(y, 1), 1) .* y, 2);

% Compute dot products between all points.
S = pdist2(V, V, fun);

% Compute shared support.
S = (acos(S) - 2*acos(h)) <= 0;

end