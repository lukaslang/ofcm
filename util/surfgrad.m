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
function f = surfgrad(Ns, c, xi)
%SURFGRAD Computes the surface gradient of a function on the sphere given 
%in spherical harmonics expansion.
%   
%   f = surfgrad(Ns, c, xi) takes list of degrees Ns, 
%   coefficients c, a radius function rho and evaluation points xi, and
%   returns the surface gradient.
%
%   Ns is a vector of consecutive non-negative integers.
%   c is a vector of length k.
%   rho is a vector of length n.
%   xi is an n-by-2 matrix of coordinates [el, az], where el in [0, pi] and
%   az in [0, 2pi).
%
%   f is a matrix of size [n, 3].
n = size(xi, 1);

% Get coordinates.
el = xi(:, 1);

% Evaluate first partial derivatives of spherical harmonics.
[d1Yij, d2Yij] = spharmderivn(Ns, xi);

% Compute partial derivatives of rho.
d1rho = d1Yij*c;
d2rho = d2Yij*c;

% Find points above certain elevation.
idx = el > pi/4;

% Compute partial derivatives of parametrisation in first chart.
[d1, d2] = sphtanbasis(xi(idx, :), eye(3));

% Compute gradient.
f1 = bsxfun(@times, d1rho(idx, :), d1) + bsxfun(@times, d2rho(idx, :) ./ (sin(el(idx)).^2), d2);

% Choose second chart.
theta = pi/2;
R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];

% Convert to spherical coordinates in new chart.
[el, az] = convertsphcoord(xi(~idx, :), R);

% Compute points on sphere.
x = sphcoord([el, az], R);

% Compute partial derivatives of parametrisation in second chart.
[d1, d2] = sphtanbasis([el, az], R);

% Compute partial derivatives via chain rule.
d1F1 = - d1(:, 3) ./ sqrt(1 - x(:, 3).^2);
d2F1 = - d2(:, 3) ./ sqrt(1 - x(:, 3).^2);
d1F2 = (d1(:, 2) .* x(:, 1) - x(:, 2) .* d1(:, 1)) ./ (x(:, 1).^2 + x(:, 2).^2);
d2F2 = (d2(:, 2) .* x(:, 1) - x(:, 2) .* d2(:, 1)) ./ (x(:, 1).^2 + x(:, 2).^2);

d1rhop = d1rho(~idx) .* d1F1 + d2rho(~idx) .* d1F2;
d2rhop = d1rho(~idx) .* d2F1 + d2rho(~idx) .* d2F2;

% Compute gradient in second map.
f2 = bsxfun(@times, d1rhop, d1) + bsxfun(@times, d2rhop ./ (sin(el).^2), d2);

% Merge gradients.
f = zeros(n, 3);
f(idx, :) = f1;
f(~idx, :) = f2;

end