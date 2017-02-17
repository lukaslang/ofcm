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
function [ofc, L] = of(Ns, cs, X, k, h, xi, w, gradf, dtdf, alpha)
%OF Computes coefficients of optical flow.
%   
%   ofc = OF(Ns, cs, x, k, h, xi, w, gradf, dtdf, alpha) takes degrees Ns
%   and Fourier coefficients cs of the surface, center points X of basis 
%   functions, degree k, a scaling factor h, quadrature rule [xi, w], and 
%   gradient gradf of data, temporal derivative dtdf, and a regularisation 
%   parameter alpha, and returns the coefficients ofc of the solution.
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
%   alpha > 0 is a scalar.
%
%   ofc is a vector of length m.
%   L is a struct with info about linear system solve.

% Set default memory 1GB.
mem = 1024^3;

% Compute optimality conditions.
[~, A, D, ~, b] = optcond(Ns, cs, X, k, h, xi, w, gradf, dtdf, ones(size(dtdf, 1), 1), mem);

% Solve linear system.
[ofc, L] = solvesystem(A + alpha * D, b, 1e-6, min(1000, size(A, 1)));

end