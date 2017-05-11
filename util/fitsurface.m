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
function [c, Y] = fitsurface(N, X, beta0, beta1, s)
%FITSURFACE Fits an evolving sphere-like surface to data.
%
%   [c, Y] = FITSURFACE(N, X, beta0, beta1) takes cell array data X in R^3 
%   and fits an evolving sphere-like surface based on spherical harmonics.
%   Non-negative scalar beta0 is a spatial and beta1 a temporal 
%   regularisation parameter and real scalar s is the parameter of the 
%   Sobolev seminorm of H^{s}(S^2, R).
%
%   X is a cell array each containing an n-by-3 matrix, N are the degrees 
%   of spherical harmonics and must be non-negative consecutive integers.
%
%   c is a cell array with vectors of coefficients of the spherical
%   harmonics Y and has length dim. Y itself is a cell array with n-by-dim 
%   matrices of scalar spherical harmonics evaluated at points X projected
%   to the unit sphere.
assert(iscell(X));
assert(isscalar(beta0));
assert(beta0 >= 0);
assert(isscalar(beta1));
assert(beta1 >= 0);
assert(isscalar(s));
assert(isvector(N));
assert(all(N >= 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

nframes = length(X);

% Compute function values (Euclidean norm of data X).
len = cellfun(@(x) sqrt(sum(x.^2, 2)), X, 'UniformOutput', false); 

% Project data to unit sphere (i.e. divide by length).
Xproj = cellfun(@(x, y) bsxfun(@rdivide, x, y), X, len, 'UniformOutput', false);

% Create diagonal matrix.
E = spharmeigs(N);
dim = length(E);
D = beta0 * diag(repmat(E .^ s, nframes, 1));

% Create matrix A and vector b.
A = sparse(nframes * dim, nframes * dim);
b = zeros(nframes * dim, 1);

Y = cell(nframes, 1);
for k=1:nframes
    % Compute spherical harmonics at projected points.
    Y{k} = spharmn(N, Xproj{k});

    % Place submatrix.
    A((k-1)*dim + 1: k*dim, (k-1)*dim + 1: k*dim) = Y{k}'*Y{k};

    % Compute vector b.
    b((k-1)*dim + 1: k*dim) = sum(bsxfun(@times, Y{k}, len{k}), 1)';
end

% Create matrix for temporal regularisation.
B = zeros(nframes*dim, nframes*dim);
if(nframes > 1)
    B = diag([beta1*ones(1, dim), 2*beta1*ones(1, (nframes-2)*dim), beta1*ones(1, dim)]);
    B = B + diag(-beta1*ones(1, (nframes-1)*dim), dim) + diag(-beta1*ones(1, (nframes-1)*dim), -dim);
end

% Solve linear system.
c = gmres(A + D + B, b, [], 1e-6, min(2000, size(A, 1)));
c = mat2cell(c, repmat(dim, nframes, 1));

end