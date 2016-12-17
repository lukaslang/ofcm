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
function [cmc, L] = cm(Ns, cs1, cs2, X, k, h, xi, w, gradf, dtdf, f, alpha)
%OFDM Computes coefficients of mass conservation.
%   
%   cmc = CM(Ns, cs1, cs2, x, k, h, xi, w, gradf, dtdf, alpha) takes degrees Ns
%   and Fourier coefficients cs of the surface, center points X of basis 
%   functions, degree k, a scaling factor h, quadrature rule [xi, w], and 
%   gradient gradf of data, temporal derivative dtdf, data f, and a 
%   regularisation parameter alpha, and returns the coefficients cmc of the
%   solution.
%
%   Ns is a vector of consecutive non-negative integers.
%   cs is a vector of Fourier coefficients (according to degrees Ns).
%   X is a matrix [m, 3] of points on the 2-sphere.
%   k is an integer degree of basis functions.
%   h in [-1, 1] is a scaling factor.
%   xi is a [n, 2] matrix of spherical coordinates [el, az] of evaluation
%   points.
%   w a quadrature weight-vector of length n.
%   gradf is a matrix of size [n, 3] of vectors.
%   dtdf is a vector of length n.
%   f is a vector of length n.
%   alpha > 0 is a scalar.
%
%   cmc is a vector of length m.
%   L is a struct with info about linear system solve.

% Compute optimality conditions.
[~, A, D, ~, b] = optcondcm(Ns, cs1, cs2, X, k, h, xi, w, gradf, dtdf, f, ones(size(f, 1), 1));

% Solve linear system.
[cmc, L] = solvesystem(A + alpha * D, b, 1e-6, 1000);

end