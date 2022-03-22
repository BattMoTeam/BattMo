classdef SolidDiffusionModel < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();

        % Physicochemical properties
        volumetricSurfaceArea  % Surface area to volume,       [m2 m^-3]
        rp                     % Particle radius               [m]
        D0                     % Diffusion coefficient         [m]
        EaD
        
        np  % Number of particles
        N   % Discretization parameters in spherical direction

    end

    methods

        function model = SolidDiffusionModel(paramobj)

            model = model@BaseModel();

             % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);

            fdnames = {'rp'                    , ...
                       'volumetricSurfaceArea' , ...
                       'EaD'                   , ...
                       'D0'                    , ...
                       'np'                    , ...
                       'N'};

            model = dispatchParams(model, paramobj, fdnames);
            model.operators = model.setupOperators();
            
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);
                        
            varnames = {};
            % concentration
            varnames{end + 1} = 'c';
            % surface concentration
            varnames{end + 1} = 'cSurface';
            % Temperature
            varnames{end + 1} = 'T';
            % Diffusion coefficient
            varnames{end + 1} = 'D';
            %
            varnames{end + 1} = 'R';
            % Mass accumulation term
            varnames{end + 1} = 'massAccum';
            % flux term
            varnames{end + 1} = 'flux';
            % Mass source term
            varnames{end + 1} = 'massSource';
            % Mass conservation equation
            varnames{end + 1} = 'massCons';
            % Mass conservation equation
            varnames{end + 1} = 'solidDiffusionEq';
            
            model = model.registerVarNames(varnames);

            fn = @SolidDiffusionModel.updateDiffusionCoefficient;
            inputnames = {'T'};
            model = model.registerPropFunction({'D', fn, inputnames});

            fn = @SolidDiffusionModel.updateFlux;
            inputnames = {'c', 'D'};
            model = model.registerPropFunction({'flux', fn, inputnames});
            
            fn = @SolidDiffusionModel.updateMassConservation;
            inputnames = {'massAccum', 'flux', 'massSource'};
            model = model.registerPropFunction({'massCons', fn, inputnames});

            fn = @SolidDiffusionModel.updateMassSource;
            model = model.registerPropFunction({'massSource', fn, {'R'}});
            
            fn = @SolidDiffusionModel.assembleSolidDiffusionEquation;
            model = model.registerPropFunction({'cSurface', fn, {'c', 'massSource', 'cSurface'}});
            
        end
        
        function operators = setupOperators(model)
            
            np = model.np;
            N  = model.N;
            rp = model.rp;
            
            celltbl.cells = (1 : np)';
            celltbl = IndexArray(celltbl);

            Scelltbl.Scells = (1 : N)';
            Scelltbl = IndexArray(Scelltbl);

            cellScelltbl = crossIndexArray(celltbl, Scelltbl, {}, 'optpureproduct', true);
            cellScelltbl = sortIndexArray(cellScelltbl, {'cells', 'Scells'});

            endScelltbl.Scells = N;
            endScelltbl = IndexArray(endScelltbl);
            endcellScelltbl = crossIndexArray(cellScelltbl, endScelltbl, {'Scells'});

            G = cartGrid(N, rp); 
            r = G.nodes.coords;

            G.cells.volumes   = 4/3*pi*(r(2 : end).^3 - r(1 : (end - 1)).^3);
            G.cells.centroids = (r(2 : end) + r(1 : (end - 1)))./2;

            G.faces.centroids = r;
            G.faces.areas     = 4*pi*r.^2;
            G.faces.normals   = G.faces.areas;

            rock.perm = ones(N, 1);
            rock.poro = ones(N, 1);

            op = setupOperatorsTPFA(G, rock);
            C = op.C;
            T = op.T;
            T_all = op.T_all;

            % We use that we know *apriori* the indexing given by cartGrid
            Tbc = T_all(N); % half-transmissibility for of the boundary face
            Tbc = repmat(Tbc, np, 1);
            
            Sfacetbl.Sfaces = (1 : (N - 1))'; % index of the internal faces (correspond to image of C')
            Sfacetbl = IndexArray(Sfacetbl);
            cellSfacetbl = crossIndexArray(celltbl, Sfacetbl, {}, 'optpureproduct', true);
            
            Grad = -diag(T)*C;

            [i, j, grad] = find(Grad);
            ScellSfacetbl.Sfaces = i;
            ScellSfacetbl.Scells = j;
            ScellSfacetbl = IndexArray(ScellSfacetbl);

            cellScellSfacetbl = crossIndexArray(celltbl, ScellSfacetbl, {}, 'optpureproduct', true);

            map = TensorMap();
            map.fromTbl = ScellSfacetbl;
            map.toTbl = cellScellSfacetbl;
            map.mergefds = {'Scells', 'Sfaces'};
            map = map.setup();

            grad = map.eval(grad);

            prod = TensorProd();
            prod.tbl1 = cellScellSfacetbl;
            prod.tbl2 = cellScelltbl;
            prod.tbl3 = cellSfacetbl;
            prod.mergefds = {'cells'};
            prod.reducefds = {'Scells'};

            Grad = SparseTensor();
            Grad = Grad.setFromTensorProd(grad, prod);
            Grad = Grad.getMatrix();

            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = cellSfacetbl;
            map.mergefds = {'cells'};
            map = map.setup();

            flux = @(D, c) - map.eval(D).*(Grad*c);

            [i, j, d] = find(C');

            clear ScellSfacetbl
            ScellSfacetbl.Scells = i;
            ScellSfacetbl.Sfaces = j;
            ScellSfacetbl = IndexArray(ScellSfacetbl);

            cellScellSfacetbl = crossIndexArray(celltbl, ScellSfacetbl, {}, 'optpureproduct', true);

            map = TensorMap();
            map.fromTbl = ScellSfacetbl;
            map.toTbl = cellScellSfacetbl;
            map.mergefds = {'Scells', 'Sfaces'};
            map = map.setup();

            d = map.eval(d);

            prod = TensorProd();
            prod.tbl1 = cellScellSfacetbl;
            prod.tbl2 = cellSfacetbl;
            prod.tbl3 = cellScelltbl;
            prod.mergefds = {'cells'};
            prod.reducefds = {'Sfaces'};
            prod = prod.setup();

            divMat = SparseTensor();
            divMat = divMat.setFromTensorProd(d, prod);
            divMat = divMat.getMatrix();

            div = @(u) divMat*u;

            %% External flux map (from the boundary conditions)

            map = TensorMap();
            map.fromTbl = endcellScelltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'cells', 'Scells'};
            map = map.setup();

            f = map.eval(ones(endcellScelltbl.num, 1));

            prod = TensorProd();
            prod.tbl1 = cellScelltbl;
            prod.tbl2 = celltbl;
            prod.tbl3 = cellScelltbl;
            prod.mergefds = {'cells'};

            mapFromBc = SparseTensor();
            mapFromBc = mapFromBc.setFromTensorProd(f, prod);
            mapFromBc = mapFromBc.getMatrix();

            mapToBc = mapFromBc';
            
            vols = G.cells.volumes;

            map = TensorMap();
            map.fromTbl = Scelltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'Scells'};
            map = map.setup();

            vols = map.eval(vols);

            operators = struct('div'      , div       , ...
                               'flux'     , flux      , ...
                               'mapFromBc', mapFromBc , ...
                               'mapToBc'  , mapToBc   , ...
                               'Tbc'      , Tbc       , ...
                               'vols'     , vols);
            
        end

        function state = updateMassSource(model, state)
            
            op = model.operators;
            rp = model.rp;
            volumetricSurfaceArea = model.volumetricSurfaceArea;
            
            R = state.R;

            R = op.mapFromBc*R;
            
            state.massSource = -R/(volumetricSurfaceArea)*(4*pi*rp^2);
            
        end
        
        function state = updateAccumTerm(model, state, state0, dt)

            op = model.operators;
            
            c = state.c;
            c0 = state0.c;
            
            state.accumTerm = 1/dt*op.vols.*(c - c0);
            
        end
        
        function state = updateMassConservation(model, state)
           
            op = model.operators;
            
            flux       = state.flux;
            massSource = state.massSource;
            accumTerm  = state.accumTerm;
            
            state.massCons = accumTerm + op.div(flux) - massSource;

        end
        
        
        
        function state = updateFlux(model, state)
            
            op = model.operators;
            
            c = state.c;
            D = state.D;
            
            state.flux = op.flux(D, c);
            
        end

        function state = updateDiffusionCoefficient(model, state)

            Tref = 298.15;  % [K]

            T = state.T;

            R   = model.constants.R;
            D0  = model.D0;
            EaD = model.EaD;

            % Calculate solid diffusion coefficient, [m^2 s^-1]
            D = D0.*exp(-EaD./R*(1./T - 1/Tref));

            state.D = D;
            
        end
    
        function state = assembleSolidDiffusionEquation(model, state)
            
            op = model.operators;

            c     = state.c;
            cSurf = state.cSurface;
            src   = state.massSource;
            
            eq = op.Tbc.*(op.mapToBc*c - cSurf) + op.mapToBc*src;
            
            state.solidDiffusionEq = eq;
            
        end
        
        
    end
    
end


%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
    

