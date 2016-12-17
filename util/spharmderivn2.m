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
function [d11, d12, d21, d22] = spharmderivn2(N, xi)
%SPHARMDERIVN2 Computes second partial derivatives of spherical harmonics.
%   
%   [d1, d2] = SPHARMDERIVN(N, xi) takes a list of degrees N >= 0, a list 
%   of evaluation points xi, and returns second partial derivatives of all 
%   spherical harmonics Y_nj of degree N.
%
%   N is a vector of non-negative consecutive integers.
%   xi is a n-by-2 matrix [el, az] where el in [0, pi] and az in [0, 2pi).
%
%   Output arguments are n-by-m matrices, where m are all orders of all 
%   degrees up to N of spherical harmonics,
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
d11 = zeros(n, dim);
d12 = zeros(n, dim);
d21 = zeros(n, dim);
d22 = zeros(n, dim);

c = 1;
for k=N
    % Generate second partial derivatives of scalar spherical harmonics.
    if(k > 0)
        [d11t, d12t, d21t, d22t] = spharmderiv2(k, xi);
        d11(:, c:(c+2*k)) = d11t;
        d12(:, c:(c+2*k)) = d12t;
        d21(:, c:(c+2*k)) = d21t;
        d22(:, c:(c+2*k)) = d22t;
    else
        d11(:, c:(c+2*k)) = zeros(n, 1);
        d12(:, c:(c+2*k)) = zeros(n, 1);
        d21(:, c:(c+2*k)) = zeros(n, 1);
        d22(:, c:(c+2*k)) = zeros(n, 1);
    end    
    c = c + 2*k + 1;
end

end