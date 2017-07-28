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
function streamlines2(P, v, S, stepsize, maxit, cmap, lineWidth)
%STREAMLINES2 Plots integral curves in the plane.
%   
%   STREAMLINES2(P, v, S, stepsize, maxit, cmap, lineWidth) takes points P
%   and a vector field v specified at P, and plots colored integral curves
%   for an artificial time. The number of steps is given by maxit.
%   S are seed points of the integral curves in R^2. cmap is used to
%   interpolate colors. lineWidth is the line width used for plotting.
%
%   Note that line segments are plotted with increasing time so that they
%   may overlap previously drawn segments.

% Initialise vertices.
verts = cell(size(S, 1), 1);

% Interpolate vector field.
F1 = scatteredInterpolant(P, v(:, 1), 'linear', 'none');
F2 = scatteredInterpolant(P, v(:, 2), 'linear', 'none');

% Integrate flow.
for k=1:length(verts)
    posx = S(k, 1);
    posy = S(k, 2);
    verts{k} = [posx, posy];
    for l=1:maxit
        x = posx + stepsize * F1(posx, posy);
        y = posy + stepsize * F2(posx, posy);
        % Stop if point has moved out of convex hull.
        if(isnan(x) || isnan(y))
            break;
        end
        % Stop if point hasen't moved.
        if(hypot(x - posx, y - posy) <= eps)
            break;
        end
        verts{k} = [verts{k}; x, y];
        posx = x;
        posy = y;
    end
end

% Interpolate colours.
c = 1:maxit;
cmaps = colormap(cmap);
y = linspace(max(c),min(c),size(cmaps,1));
cm = spline(y, cmaps', c);

% Plot lines.
for l=1:maxit-1
    X = [];
    Y = [];
    for k=1:length(verts)
        v = verts{k};
        if(l >= size(v, 1))
            continue;
        end        
        X = [X; v(l, 1), v(l+1, 1)];
        Y = [Y; v(l, 2), v(l+1, 2)];
    end
    if(isempty(X) || isempty(Y))
        break;
    end
    line(X', Y', 'color', cm(:, l), 'LineWidth', lineWidth);
end
end