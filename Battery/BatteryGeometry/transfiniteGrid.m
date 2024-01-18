function G = transfiniteGrid(u, v, f1, f2, g1, g2, varargin)
% Construct a 2D grid using transfinite interpolation
%
% SYNOPSIS:
%   G = transfiniteGrid(u, v, f1, f2, g1, g2)
%
% DESCRIPTION:
%
%   Create a grid using transfinite interpolation of boundary curves.
%   In 2D there are 4 curves: f1(u) and f2(u) where u \in [0,1], as
%   well as g1(v) and g2(v) where v \in [0,1]. The topology is the
%   same as for cartGrid.
%
%   We allow for constructing grids with f1==f2 by setting the flag
%   "closed" to "true". If a hole is formed, this can be filled with a
%   single cell using the flag "fill" set to "true" (setting "fill" to
%   "true" implies setting "closed" to "true").
%
%   Depending on the curves f1, f2, g1, g2 the grid may be arbitrarily
%   twisted. No checks for this are made.
%
%          f2
%        -----
%        |   |
%     g1 |   | g2
%        |   |
%        -----
%          f1
%
% PARAMETERS:
%   u        - parameter for f1 and f2. u \in [0,1]
%   v        - parameter for g1 and g2. v \in [0,1]
%   f1, f2   - curves opposite each other parameterized by u
%   g1, g2   - curves opposite each other parameterized by v
%
% OPTIONAL PARAMETERS:
%   blending - Set the blending function order (only linear is
%              supported)
%   closed   - Create a closed grid if f1==f2
%   fill     - Fill a possible hole in case of a closed grid (implies
%              closed=true)
%   plot     - Plot node, face and cell numbers for debugging
%
% RETURNS:
%   G - An MRST grid
%
% EXAMPLE:
%  % Create a grid of a circle quadrant
%  r = 10;
%  R = 20;
%  Ns = 5;
%  s = linspace(0, 1, Ns)';
%  f1 = [(1-s)*r + s*R, zeros(size(s))];
%  f2 = [zeros(size(s)), (1-s)*r + s*R];
%  Nt = 7;
%  t = linspace(0, 1, Nt)';
%  g1 = r*[cos(t*pi/2), sin(t*pi/2)];
%  g2 = R*[cos(t*pi/2), sin(t*pi/2)];
%  G = transfiniteGrid(s, t, f1, f2, g1, g2);
%  figure, plotGrid(G), axis equal
%  text(mean(f1(:,1)), mean(f1(:,2)), 'f1')
%  text(mean(f2(:,1)), mean(f2(:,2)), 'f2')
%  text(mean(g1(:,1)), mean(g1(:,2)), 'g1')
%  text(median(g2(:,1)), mean(g2(:,2)), 'g2')
%
%
%  % Create a grid of a circle
%  r = 1;
%  R = 20;
%  Ns = 10;
%  s = linspace(0, 1, Ns)';
%  f1 = [(1-s)*r + s*R, zeros(size(s))];
%  f2 = f1;
%  Nt = 50;
%  t = linspace(0, 1, Nt)';
%  g1 = r*[cos(t*2*pi), sin(t*2*pi)];
%  g2 = R*[cos(t*2*pi), sin(t*2*pi)];
%  G = transfiniteGrid(s, t, f1, f2, g1, g2, 'fill', true);
%  figure, plotGrid(G), axis equal
%
% SEE ALSO:
%   `cartGrid`, `tensorGrid`, `grid_structure`, `computeGeometry`

%{
  Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

  This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

  MRST is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  MRST is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('blending', 'linear', ...
                 'closed'  , false, ...
                 'fill'    , false, ...
                 'plot'    , false);
    opt = merge_options(opt, varargin{:});

    if opt.fill
        opt.closed = true;
    end

    if any(u<0 | u>1)
        warning('u should be in [0,1]')
    end
    if any(v<0 | v>1)
        warning('v should be in [0,1]')
    end

    if contains(opt.blending, 'linear')
        a0 = 1-u;
        a1 = u;
        b0 = 1-v;
        b1 = v;
    else
        error('Unknown blending %s; only linear blending is implemented', opt.blending);
    end

    Nu = numel(u);
    Nv = numel(v);
    coords = zeros(Nu*Nv, 2);
    cnt = 1;

    % To be vectorized
    for iv = 1:Nv
        for iu = 1:Nu
            c1 = a0(iu)*g1(iv,:) + a1(iu)*g2(iv,:);
            c2 = b0(iv)*f1(iu,:) + b1(iv)*f2(iu,:);
            c3 = a0(iu)*b0(iv)*g1(1,:) + a0(iu)*b1(iv)*g1(end,:) + a1(iu)*b0(iv)*f1(end,:) + a1(iu)*b1(iv)*g2(end,:);
            coords(cnt,:) = c1 + c2 - c3;
            cnt = cnt+1;
        end
    end

    % Grid topology from cartGrid
    G = cartGrid([Nu-1, Nv-1]);
    G.cells.faces(:,2) = [];
    G.faces = rmfield(G.faces, 'tag');
    G.type = {'transfinite'};
    G.nodes.coords = coords;
    findNeighbors = false;

    if opt.closed

        % If we have a closed curve where f1 == f2, not g1 == g2
        assert(norm(f1-f2, 'inf') < norm(g1-g2, 'inf'), 'only f1==f2 is implemented');
        findNeighbors = true;
        G.faces = rmfield(G.faces, 'neighbors');

        % Remove last Nu nodes
        nodes = G.nodes.num + (-Nu+1:0);
        G.nodes.coords(nodes, :) = [];
        G.nodes.num = G.nodes.num - numel(nodes);

        % Map nodes to 1:Nu
        for k = 1:Nu
            idx = G.faces.nodes == nodes(k);
            G.faces.nodes(idx) = k;
        end

        % Remove the last edges
        edges = G.faces.num + (-numel(nodes)+2:0);
        G.faces.nodes(end-2*numel(edges)+1:end) = [];
        G.faces.num = G.faces.num - numel(edges);
        G.faces.nodePos = G.faces.nodePos(1:G.faces.num+1);

        % Map edges. Since f1==f2 the new edge numbers are
        % Nu*(Nv-1)+(1:Nu-1)
        for k = 1:numel(edges)
            idx = G.cells.faces == edges(k);
            G.cells.faces(idx) = Nu*(Nv-1) + k;
        end

    end

    if opt.fill

        % Fill the hole with a single cell
        assert(opt.closed)
        G.getNumberOfCells() = G.getNumberOfCells() + 1;
        f = 1:Nu:(Nv-1)*Nu;
        G.cells.faces = [G.cells.faces; f'];
        G.cells.facePos(end+1) = G.cells.facePos(end) + numel(f);

    end

    G = computeGeometry(G, 'findNeighbors', findNeighbors);

    if opt.plot
        debugPlot(G);
    end

end

function debugPlot(G)

    figure, hold on, axis equal tight, grid on
    plotGrid(G, 'facealpha', 0.5)
    plot(G.nodes.coords(:,1), G.nodes.coords(:,2), 'k.');
    text(G.nodes.coords(:,1), G.nodes.coords(:,2), num2str((1:G.nodes.num)'), 'color', 'k');
    plot(G.faces.centroids(:,1), G.faces.centroids(:,2), 'r.')
    text(G.faces.centroids(:,1), G.faces.centroids(:,2), num2str((1:G.faces.num)'), 'color', 'r');
    plot(G.cells.centroids(:,1), G.cells.centroids(:,2), 'b.')
    text(G.cells.centroids(:,1), G.cells.centroids(:,2), num2str((1:G.getNumberOfCells())'), 'color', 'b');

end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
