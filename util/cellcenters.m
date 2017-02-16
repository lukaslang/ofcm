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
function C = cellcenters(f, sigma, size, t, area)
%CELLCENTERS Returns approximate cell centres.
%   
%   C = CELLCENTERS(f, sigma, size, t, area) takes a cell array of
%   volumetric data and returns approximate cell centers characterized by
%   local maxima in intensity. t and area are thresholding parameters for
%   intensity and area of the local maximum, respectively.

P = cell(length(f), 1);
for k=1:length(f)
    % Apply Gaussian filter to frame.
    fg = imgaussian(double(f{k}), sigma, size);
    % Compute regionprops of local maxima.
    P{k} = regionprops(imregionalmax(fg, 26), fg, 'Area', 'Centroid', 'MaxIntensity');
end

% Select regions by thresholding.
I = cellfun(@(p) [p.MaxIntensity] > t & [p.Area] <= area, P, 'UniformOutput', false);

% Apply selection.
C = cellfun(@(p, idx) reshape([p(idx).Centroid], 3, numel(find(idx)))', P, I, 'UniformOutput', false);

end