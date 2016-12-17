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
function [u, L] = solvesystem(A, b, tol, maxit)
%SOLVESYSTEM Solves a linear system.
%
%   [u, L] = SOLVESYSTEM(A, b, maxit) solves the linear system 
%   A*u = b. tol is the desired relative residual. maxit is the maximum 
%   number of iterations during the linear system solve.
%
%   A is a matrix of size [n, n].
%   b a vector of length n.
%
%   u is a vector of length n.
%
%   L is a struct containing information about the linear system solve.

assert(isscalar(maxit));
assert(maxit > 0);

% Store norm of rhs.
L.rhs = norm(b, 2);

% Solve linear system.
ticId = tic;
[u, flag, relres, iter, resvec] = gmres(A, b, [], tol, maxit);

% Store solver information.
L.time = toc(ticId);
L.flag = flag;
L.relres = relres;
L.iter = iter;
L.resvec = resvec;
L.tol = tol;
L.maxit = maxit;
L.solver = 'gmres';
L.restart = 0;

end