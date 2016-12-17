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
function [d1, d2] = spharmderivn(N, xi)
%SPHARMDERIVN Computes partial derivatives of spherical harmonics.
%   
%   [d1, d2] = SPHARMDERIVN(N, xi) takes a list of degrees N >= 0, a list 
%   of evaluation points xi, and returns partial derivatives of all 
%   spherical harmonics Y_nj of degree N.
%
%   N is a vector of non-negative consecutive integers.
%   xi is a n-by-2 matrix [el, az] where el in [0, pi] and az in [0, 2pi).
%
%   d1 and d2 are n-by-m matrices, where m are all orders of all degrees up
%   to N of spherical harmonics,
%   i.e. m = (N(end)^2 + 2*N(end) - N(1)^2 + 1).

% Check if N is an interval of consecutive non-negative integers.
assert(isvector(N));
assert(all(N >= 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

n = size(xi, 1);

% Compute the dimension.
dim = (N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Create scalar spherical harmonics.
d1 = zeros(n, dim);
d2 = zeros(n, dim);

c = 1;
for k=N
    % Generate partial derivatives of scalar spherical harmonics of degree k.
    [d1t, d2t] = spharmderiv(k, xi);
    d1(:, c:(c+2*k)) = d1t;
    d2(:, c:(c+2*k)) = d2t;  
    c = c + 2*k + 1;
end

end