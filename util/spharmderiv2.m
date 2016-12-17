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
function [d11, d12, d21, d22] = spharmderiv2(N, xi)
%SPHARMDERIV2 Computes second partial derivatives of spherical harmonics.
%   
%   [d11, d12, d21, d22] = SPHARMDERIV2(N, xi) takes a degree N, a list of 
%   evaluation points xi, and returns partial derivatives of all spherical
%   harmonics of degree N.
%
%   N is a positive integer.
%   xi is a n-by-2 matrix [el, az] where el in [0, pi] and az in [0, 2pi).
%
%   All output arguments are n-by-m matrices, where m are all orders of 
%   degree N of spherical harmonics.
n = size(xi, 1);

% Get coordinates.
el = xi(:, 1);
az = xi(:, 2);

% Compute legendre polynomials.
P = legendre(N, cos(el), 'norm')';

% Compute polynomials.
ord = abs(-N:N);
factor1 = 0.25 * sqrt((N + ord) .* (N - ord + 1) .* (N + ord - 1) .* (N - ord + 2)) .* [repmat(2, 1, N), sqrt(2), repmat(2, 1, N)] ./ sqrt(4*pi);
factor2 = 0.5 * (ord.^2 - N^2 - N) .* [repmat(2, 1, N), sqrt(2), repmat(2, 1, N)] ./ sqrt(4*pi);
factor3 = 0.25 * sqrt((N + ord + 1) .* (N + ord + 2) .* (N - ord - 1) .* (N - ord)) .* [repmat(2, 1, N), sqrt(2), repmat(2, 1, N)] ./ sqrt(4*pi);
P1 = [zeros(n, 2), P(:, 1:end-2)];
P3 = [P(:, 3:end), zeros(n, 2)];
Np = bsxfun(@times, factor1, [flip(P1(:, 2:end), 2), P1]) + bsxfun(@times, factor2, [flip(P(:, 2:end), 2), P]) + bsxfun(@times, factor3, [flip(P3(:, 2:end), 2), P3]);

% Reorder so that order is -n,...,n.
C = repmat(-N:N, n, 1) .* repmat(az, 1, 2*N + 1);
C = [cos(C(:, 1:N + 1)), sin(C(:, N + 2:2*N + 1))];
d11 = C .* Np;

% Compute spherical harmonics.
Y = spharmp(N, az, cos(el));
d22 = - bsxfun(@times, (-N:N).^2, Y);

% Compute polynomials.
factor1 = 0.5 * sqrt((N + ord) .* (N - ord + 1)) .* [repmat(2, 1, N), sqrt(2), repmat(2, 1, N)] ./ sqrt(4*pi);
factor2 = 0.5 * sqrt((N + ord) .* (N - ord)) .* [repmat(2, 1, N), sqrt(2), repmat(2, 1, N)] ./ sqrt(4*pi);
P1 = [zeros(n, 1), P(:, 1:end-1)];
P2 = [P(:, 2:end), zeros(n, 1)];
Np = bsxfun(@times, factor1, [flip(P1(:, 2:end), 2), P1]) - bsxfun(@times, factor2, [flip(P2(:, 2:end), 2), P2]);

% Reorder so that order is -n,...,n.
C = repmat(-N:N, n, 1) .* repmat(az, 1, 2*N + 1);
C = [-sin(C(:, 1:N + 1)), cos(C(:, N + 2:2*N + 1))];
d12 = bsxfun(@times, C .* Np, -N:N);

d21 = d12;

end