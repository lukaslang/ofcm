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
function [dim, A, D, E, G, b] = optcondcm(Ns, cs1, cs2, X, k, h, xi, w, gradf, dtdf, f, s, mem)
%OPTCONDCM Computes the optimality conditions and returns a linear 
%system used in mass conservation on a sphere-like surface.
%
%   [dim, A, D, E, G, b] = OPTCONDCM(Ns, cs1, cs2, X, k, h, xi, w, gradf, dtdf, f, s) 
%   takes degrees Ns and Fourier coefficients cs1, cs2 of the surface, 
%   center points X of basis functions, degree k, a scaling factor h, 
%   quadrature rule [xi, w], and  gradient gradf of data, temporal 
%   derivative dtdf, data f, and a segmentation s, and returns the linear 
%   system representing the optimality conditions.
%
%   [dim, A, D, E, G, b] = OPTCONDCM(Ns, cs1, cs2, X, k, h, xi, w, gradf, dtdf, f, s, mem)
%   uses mem > 0 bytes for matrix multiplication.
%
%   Ns is a vector of consecutive non-negative integers.
%   cs1, cs2 are vectors of Fourier coefficients (according to degrees Ns).
%   X is a matrix [m, 3] of points on the 2-sphere.
%   k is an integer degree of basis functions.
%   h in [-1, 1] is a scaling factor.
%   xi is a [n, 2] matrix of spherical coordinates [el, az] of evaluation
%   points.
%   w a quadrature weight-vector of length n.
%   gradf is a matrix of size [n, 3] of vectors.
%   dtdf is a vector of length n.
%   f is a vector of length n.
%
%   Note that degrees Ns must be a vector of consecutive positive integers!
%   
%   The linear system returned is
%
%   (A + alpha * D + beta * E + gamma * G) * x = b.
%
%   A, D, E, and G are matrices of size [2*m, 2*m].
%   b is the right hand side and is of length 2*m.

% Number of basis functions.
dim = 2*size(X, 1);

% Number of evaluation points.
n = size(xi, 1);

% Check for memory option.
if(nargin == 13)
    assert(mem > 0);
    mem1 = mem;
    mem2 = mem;
else
    mem1 = 8*n*dim/2;
    mem2 = 3*8*n*dim/2;
end

% Compute evaluation points.
Y = sphcoord(xi, eye(3));

% Compute synthesis of evaluation points.
[~, rho] = surfsynth(Ns, Y, cs2);

% Evaluate basis functions.
[bfc1, bfc2] = vbasiscompmem(k, h, X, Y, mem1);

% Compute coordinate basis at evaluation points.
[d1, d2] = surftanbasis(Ns, cs2, xi);

% Compute derivatives of coordinates basis at evaluation points.
[d11, d12, ~, d22] = surftanbasisderiv(Ns, cs2, xi);

% Compute surface divergence of basis functions.
sdiv = surfdivmem(Ns, cs2, xi, k, h, X, mem);

% Compute sum of inner products with gradient of data and data times div.
ip = sparse(bsxfun(@times, bfc1, dot(d1, gradf, 2)') + bsxfun(@times, bfc2, dot(d2, gradf, 2)') + bsxfun(@times, sdiv, f'));

% Evaluate partial derivatives of basis functions.
[bfcd{1, 1}, bfcd{1, 2}, bfcd{2, 1}, bfcd{2, 2}] = vbasiscompderivmem(k, h, X, Y, mem2);

% Evaluate norm squared of surface gradient of rho.
gradrhosquared = surfgradnormsquared(Ns, cs2, xi);

% Evaluate metric.
[g{1,1}, g{1,2}, g{2,1}, g{2,2}] = surfmetric(Ns, cs2, rho, xi);

% Compute determinant of metric.
detg = g{1,1} .* g{2,2} - g{1,2} .* g{2,1};

% Evaluate inverse of metric.
[ginv{1,1}, ginv{1,2}, ginv{2,1}, ginv{2,2}] = surfinvmetric(Ns, cs2, rho, xi, detg);

% Compute Christoffel symbols with respect to space.
[c{1,1,1}, c{1,1,2}, c{1,2,1}, c{1,2,2}, c{2,2,1}, c{2,2,2}] = surfchristoffel(Ns, cs2, rho, xi);
c{2, 1, 1} = c{1, 2, 1};
c{2, 1, 2} = c{1, 2, 2};

% Compute scaling factor for manifold integration.
intf = w .* rho .* sqrt(gradrhosquared + rho.^2);

% Compute surface normals.
N = surfnormals(Ns, cs2, xi);

% Compute scalar normal part of surface velocity.
[~, dtrho] = surfsynth(Ns, Y, cs2 - cs1);
V = dtrho .* dot(Y, N, 2);

% Compute dot product of surface velocity with gradient of data.
ipvtan = dot(bsxfun(@times, Y, dtrho), gradf, 2);

% Compute total curvature.
K = (dot(d11, N, 2) .* g{2, 2} + dot(d22, N, 2) .* g{1, 1} - 2 * dot(d12, N, 2) .* g{1, 2}) ./ detg;

% Create linear system.
A = ip*spdiags(intf, 0, n, n)*(ip');

% Compute components D_i^j of covariant derivatives of basis functions.
C = cell(2, 2);
for i=1:2
    for k=1:2
        C{i, k} = sparse(bfcd{k, i} + bsxfun(@times, bfc1, c{i, 1, k}') + bsxfun(@times, bfc2, c{i, 2, k}'));
    end
end

% Compute regularisation matrix.
D = sparse(dim, dim);
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                D = D + C{i, k} * spdiags(g{k, l} .* ginv{i, j} .* intf .* s, 0, n, n) * (C{j, l}');
            end
        end
    end
end

% Create second regularisation matrix.
E = bfc1*spdiags(g{1, 1} .* intf .* (1 - s), 0, n, n)*(bfc1');
E = E + bfc1*spdiags(g{1, 2} .* intf .* (1 - s), 0, n, n)*(bfc2');
E = E + bfc2*spdiags(g{2, 1} .* intf .* (1 - s), 0, n, n)*(bfc1');
E = E + bfc2*spdiags(g{2, 2} .* intf .* (1 - s), 0, n, n)*(bfc2');

% Create divergence penalty matrix.
G = sdiv*spdiags(intf .* (1 - s), 0, n, n)*sdiv';

% Create right-hand side.
b = -ip*spdiags(intf, 0, n, n)*(dtdf - ipvtan - f .* K .* V);

end