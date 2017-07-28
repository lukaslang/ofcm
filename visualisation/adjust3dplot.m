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
function adjust3dplot
%ADJUST3DPLOT Sets axes, etc. of a 3D plot.

daspect([1, 1, 1]);
set(gca, 'ZLim', [0, 420]);
set(gca, 'XLim', [-400, 400]);
set(gca, 'YLim', [-400, 400]);
set(gca, 'XTick', -450:150:450);
set(gca, 'YTick', -450:150:450);
set(gca, 'ZTick', -450:150:450);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 14);

end