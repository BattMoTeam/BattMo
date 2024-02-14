classdef GenericGrid
%
% Grid structure that can be used for a parent grid in a sub-grid (see SubGrid class).
%
% The geometrical properties are all stored in a immutable handle (see TwoPointFiniteVolumeGeometry class).
%
% By using handle, we can share the grid as a parent grid between different sub-grid instances.
%
% By using immutable properties, we avoid any collateral effects arising from the use of handler.
%
% The assembly is done using the function listed here that are at the end. Note that, typically, we will used the
% sub-grid versions of those as this class is mainly meant to be used as parent grid
%
% - mrstFormat
% - getVolumes
% - getNumberOfCells
% - getNumberOfFaces
% - getIntFaces
% - getFaceAreas
% - getGrad
% - getDiv
% - getTransHarmFace
% - getBcTrans
% - getTransBcHarmFace
% - getTrans
% - getFaceUpstreamValue
% - getCellFluxNorm




    properties (SetAccess = protected)

        % Description of topology using MRST grid structures: The fieds are
        % - topology.cells.facePos
        % - topology.cells.faces
        % - topology.cells.num
        % - topology.faces.nodePos
        % - topology.faces.nodes
        % - topology.faces.num
        % - topology.faces.neighbors
        % - topology.nodes.num
        % - topology.griddim
        topology

        % Helper structures that are used in the computation of the gradient, divergence, harmonic face averages, ...
        % - helpers.intfaces              (Index of the internal faces)
        % - helpers.diffop.grad           (sparse matrix used in getGradient)
        % - helpers.diffop.div            (sparse matrix used in getDiv)
        % - helpers.trans.D               (sparse matrix used in getTransHarmFace method)
        % - helpers.trans.P               (sparse matrix used in getTransHarmFace method)
        % - helpers.trans.S               (sparse matrix used in getTransHarmFace method)
        % - helpers.extfaces.faces        (index of the external faces)
        % - helpers.extfaces.cells        (index of the corresponding cells)
        % - helpers.extfaces.cellfacemap  (index of the corresponding cell-face index. This is used to retrieve half transmissibilities)
        % - helpers.faceextfacemap        (mapping from face to extface)
        helpers

        % Operators to compute norm of the flux velocity at the cell centers, see getCellFluxNorm
        % - cellFluxOperators.P
        % - cellFluxOperators.S
        cellFluxOperators


    end


    methods

        function grid = GenericGrid(G)
        % We initialize the grid using a MRST grid structure

            topology.cells.facePos   = G.cells.facePos;
            topology.cells.faces     = G.cells.faces;
            topology.cells.num       = G.cells.num;
            topology.faces.nodePos   = G.faces.nodePos;
            topology.faces.nodes     = G.faces.nodes;
            topology.faces.num       = G.faces.num;
            topology.faces.neighbors = G.faces.neighbors;
            topology.nodes.num       = G.nodes.num;
            topology.griddim         = G.griddim;

            grid.topology = topology;
            
        end


        function grid = setupCellFluxOperators(grid)
        % Helper operators for CellFluxNorm. Note : those have no AD compliant version for the moment.
            G = grid.mrstFormat();
            grid.cellFluxOperators = getCellFluxOperatorsAll(G);

        end

        function tpfvGeometry = getTPFVgeometry()
            error('abstract method')
        end


        function G = mrstFormat(grid)
        % Returns MRST grid that can be used for plotting

            G  = grid.topology;
            tg = grid.tPFVgeometry;

            d = G.griddim;

            G.type = 'generic';

            G.cells.centroids = reshape(tg.cells.centroids, d, [])';
            G.cells.volumes   = tg.cells.volumes;
            G.faces.centroids = reshape(tg.faces.centroids, d, [])';
            G.faces.normals   = reshape(tg.faces.normals, d, [])';
            G.faces.areas     = tg.faces.areas;
            G.nodes.coords    = reshape(tg.nodes.coords, d, [])';

        end


        function vols = getVolumes(grid)

            tpfaGeometry = grid.getTPFVgeometry();
            vols = tpfaGeometry.cells.volumes;

        end

        function nc = getNumberOfCells(grid)

            nc = grid.topology.cells.num;

        end


        function nf = getNumberOfFaces(grid)

            nf = grid.topology.faces.num;

        end

        function intfaces = getIntFaces(grid)

            intfaces = grid.helpers.intfaces;

        end

        function areas = getFaceAreas(grid)

            tpfaGeometry = grid.getTPFVgeometry();
            areas = tpfaGeometry.faces.areas;

        end

        function v = getGrad(grid, c)

            v = grid.helpers.diffop.grad*c;

        end

        function u = getDiv(grid, v)

            u = grid.helpers.diffop.div*v;

        end

        function v = getFaceUpstreamValue(grid, flag, x)

            intfaces = grid.helpers.intfaces;
            N        = grid.topology.faces.neighbors;
            nc       = grid.getNumberOfCells();
            nf       = numel(intfaces);

            v = faceUpstr(flag, x, N(intfaces, :), [nf, nc]);

        end

        function u = getTransHarmFace(grid, c)
        % Returns fluxes for each internal faces for the cell-valued vector c

            op           = grid.helpers.trans;
            tpfvGeometry = grid.getTPFVgeometry();
            hT           = tpfvGeometry.hT;

            u = 1 ./ (op.S * ( 1 ./ ((op.D*c) .* (op.P*hT))));

        end

        function [bchT, bccells, bcsgn] = getBcTrans(grid, bcfaces)

            tpfvGeometry = grid.getTPFVgeometry();
            hT           = tpfvGeometry.hT;

            exf = grid.helpers.extfaces;

            extfaceind = grid.helpers.faceextfacemap(bcfaces);

            bccells = exf.cells(extfaceind);
            bcsgn   = exf.sgn(extfaceind);
            hTind   = exf.cellfacemap(extfaceind);

            bchT = hT(hTind);

        end

        function [bchT, bccells, bcsgn] = getTransBcHarmFace(grid, c, bcfaces)

            [bchT, bccells, bcsgn] = grid.getBcTrans(bcfaces);

            bchT = bchT.*c(bccells);

        end

        function T = getTrans(grid)

            tpfvGeometry = grid.getTPFVgeometry();
            T = tpfvGeometry.T;

        end


        function jsq = getCellFluxNorm(grid, u)

            P = grid.cellFluxOperators.P;
            S = grid.cellFluxOperators.S;

            j = P*u;
            jsq = j.^2;
            jsq = S*jsq;

        end


    end

    methods(Static)

        function value = applyProduct(matrixOperator, arg1, arg2)
        % Function that assembly product structure (two dispatches followed by one reduction)

            D1 = matrixOperator.D1;
            D2 = matrixOperator.D2;
            S  = matrixOperator.S;

            value = S*((D1*arg1).*(D2*arg2));

        end
    end

end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
