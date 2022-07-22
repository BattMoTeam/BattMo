classdef SolidElectrodeInterface < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();
        
        molecularWeight % SEI molecular weight [kg/mol]
        density         % SEI densisity [kg/m^3]
        D               % SEI diffusion coefficient [m^2/s]
        
        np % number of particles
        N  % discretization parameter for SEI model (planar assumption)
        
    end

    methods

        function model = SolidElectrodeInterface(paramobj)

            model = model@BaseModel();

             % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);

            fdnames = {'molecularWeight',
                       'density',
                       'D', 
                       'np',
                       'N'};
            model = dispatchParams(model, paramobj, fdnames);
            
            model.operators = model.setupOperators();
            
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);
                        
            varnames = {};
            % Solvent concentration in the SEI film
            varnames{end + 1} = 'c';
            % solvent concentration at inteface
            varnames{end + 1} = 'cInterface';
            % external surface concentration
            varnames{end + 1} = 'cExternal';
            % surface concentration
            %% FIXME : change name delta to more explicit name (sei width)
            varnames{end + 1} = 'delta';
            % SEI growth velocity
            varnames{end + 1} = 'v';
            % Reaction rate
            varnames{end + 1} = 'R';
            % Mass accumulation term
            varnames{end + 1} = 'massAccum';
            % flux term
            varnames{end + 1} = 'flux';
            % Mass source term
            varnames{end + 1} = 'massSource';
            % Mass conservation equation
            varnames{end + 1} = 'massCons';
            % flux-concentration equation at interface boundary
            varnames{end + 1} = 'interfaceBoundaryEq';
            % SEI Width equation
            varnames{end + 1} = 'widthEq';
            
            model = model.registerVarNames(varnames);

            fn = @SolidElectrodeInterface.updateFlux;
            inputnames = {'c', 'v', 'delta'};
            model = model.registerPropFunction({'flux', fn, inputnames});
            
            fn = @SolidElectrodeInterface.updateSEIgrowthVelocity;
            inputnames = {'R'};
            model = model.registerPropFunction({'v', fn, inputnames});
            
            fn = @SolidElectrodeInterface.updateMassConservation;
            inputnames = {'massAccum', 'flux', 'massSource'};
            model = model.registerPropFunction({'massCons', fn, inputnames});

            fn = @SolidElectrodeInterface.updateMassSource;
            model = model.registerPropFunction({'massSource', fn, {'R', 'c', 'v', 'cExternal'}});
            
            fn = @SolidElectrodeInterface.assembleInterfaceBoundaryEquation;
            model = model.registerPropFunction({'interfaceBoundaryEq', fn, {'c', 'cInterface', 'massSource'}});
            
            fn = @SolidElectrodeInterface.assembleWidthEquation;
            model = model.registerPropFunction({'widthEq', fn, {'delta', 'v'}});

            fn = @SolidElectrodeInterface.updateMassAccumTerm;
            model = model.registerPropFunction({'massAccum', fn, {'delta', 'c'}});

        end
        
        function operators = setupOperators(model)
            
            np = model.np;
            N  = model.N;
            D  = model.D;
            
            celltbl.cells = (1 : np)';
            celltbl = IndexArray(celltbl);

            Scelltbl.Scells = (1 : N)';
            Scelltbl = IndexArray(Scelltbl);

            cellScelltbl = crossIndexArray(celltbl, Scelltbl, {}, 'optpureproduct', true);
            cellScelltbl = sortIndexArray(cellScelltbl, {'cells', 'Scells'});

            endScelltbl.Scells = N;
            endScelltbl = IndexArray(endScelltbl);
            endcellScelltbl = crossIndexArray(cellScelltbl, endScelltbl, {'Scells'});
            
            startScelltbl.Scells = 1;
            startScelltbl = IndexArray(startScelltbl);
            startcellScelltbl = crossIndexArray(cellScelltbl, startScelltbl, {'Scells'});
            
            Sfacetbl.Sfaces = (1 : (N - 1))'; % index of the internal faces (correspond to image of C')
            Sfacetbl = IndexArray(Sfacetbl);
            cellSfacetbl = crossIndexArray(celltbl, Sfacetbl, {}, 'optpureproduct', true);
            
            % The governing equation is setup in a fixed mesh of length 1            
            G = cartGrid(N, 1); 
            G = computeGeometry(G);
            
            rock.perm = ones(N, 1);
            rock.poro = ones(N, 1);
            op = setupOperatorsTPFA(G, rock);
            C = op.C;
            T = op.T; % in Sfacetbl
            T_all = op.T_all;
            
            TExtBc = D*T_all(N); % half-transmissibility for of the external boundary face
            TIntBc = D*T_all(1); % half-transmissibility for of the internal boundary face (close to solid particle)
            
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

            grad = map.eval(grad); % in cellScellSfacetbl

            prod = TensorProd();
            prod.tbl1 = cellScellSfacetbl;
            prod.tbl2 = cellScelltbl;
            prod.tbl3 = cellSfacetbl;
            prod.mergefds = {'cells'};
            prod.reducefds = {'Scells'};

            Grad = SparseTensor();
            Grad = Grad.setFromTensorProd(grad, prod);
            Grad = Grad.getMatrix();

            diffFlux = @(c) - D.*(Grad*c);

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

            mapFromExtBc = SparseTensor();
            mapFromExtBc = mapFromExtBc.setFromTensorProd(f, prod);
            mapFromExtBc = mapFromExtBc.getMatrix();

            mapToExtBc = mapFromExtBc';
            
            % same procedure for internal bc
            map = TensorMap();
            map.fromTbl = startcellScelltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'cells', 'Scells'};
            map = map.setup();

            f = map.eval(ones(startcellScelltbl.num, 1));

            prod = TensorProd();
            prod.tbl1 = cellScelltbl;
            prod.tbl2 = celltbl;
            prod.tbl3 = cellScelltbl;
            prod.mergefds = {'cells'};

            mapFromIntBc = SparseTensor();
            mapFromIntBc = mapFromIntBc.setFromTensorProd(f, prod);
            mapFromIntBc = mapFromIntBc.getMatrix();

            mapToIntBc = mapFromIntBc';

            %% setup Average mapping
            
            M = op.M; 
            
            [i, j, M]  = find(M);

            clear ScellSfacetbl
            ScellSfacetbl.Scells = j;
            ScellSfacetbl.Sfaces = i;
            ScellSfacetbl = IndexArray(ScellSfacetbl);

            cellScellSfacetbl = crossIndexArray(celltbl, ScellSfacetbl, {}, 'optpureproduct', true);

            map = TensorMap();
            map.fromTbl = ScellSfacetbl;
            map.toTbl = cellScellSfacetbl;
            map.mergefds = {'Scells', 'Sfaces'};
            map = map.setup();

            M = map.eval(M);
            
            prod = TensorProd();
            prod.tbl1 = cellScellSfacetbl;
            prod.tbl2 = cellScelltbl;
            prod.tbl3 = cellSfacetbl;
            prod.reducefds = {'Scells'};
            prod.mergefds = {'cells'};
            prod = prod.setup();
            
            faceAverage = SparseTensor();
            faceAverage = faceAverage.setFromTensorProd(M, prod);
            faceAverage = faceAverage.getMatrix();
            
            
            %% mapToSei (mapping from celltbl to cellScelltbl) and Centroids values
            
            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'cells'};
            map = map.setup();

            mapToSei = SparseTensor();
            mapToSei = mapToSei.setFromTensorMap(map);
            mapToSei = mapToSei.getMatrix();
            
            %%  assemble cell centroids and volumes
            
            map = TensorMap();
            map.fromTbl = Scelltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'Scells'};
            map = map.setup();
            
            cellCentroids = map.eval(G.cells.centroids);
            vols = map.eval(G.cells.volumes);
            
            operators = struct('div'          , div          , ...
                               'diffFlux'     , diffFlux     , ...
                               'faceAverage'  , faceAverage  , ...
                               'mapFromExtBc' , mapFromExtBc , ...
                               'mapToExtBc'   , mapToExtBc   , ...
                               'mapFromIntBc' , mapFromIntBc , ...
                               'mapToIntBc'   , mapToIntBc   , ...
                               'TIntBc'       , TIntBc       , ...
                               'TExtBc'       , TExtBc       , ...
                               'mapToSei'     , mapToSei     , ...
                               'cellCentroids', cellCentroids, ...
                               'vols'         , vols);
            
        end

        function state = updateMassSource(model, state)
            
            op = model.operators;
            
            R     = state.R;
            c     = state.c;
            cext  = state.cExternal;
            delta = state.delta;
            
            c = op.mapToExtBc*c;
            %% FIXME : add epsilon term before cext (as in Safari paper)
            % The convection term vanishes for xi = 1
            srcExternal = -op.TExtBc.*(c - cext);
            srcExternal = op.mapFromExtBc*srcExternal;
            
            % Here, we use conversion relation between source terms given in moving and fixed mesh : 
            % xi*v*delta*c_fixed = delta*src_moving - src_fixed (see doc).
            srcInterface = op.mapFromIntBc*(delta.*R);
            state.massSource = srcExternal + srcInterface;
            
        end

        
        function state = assembleWidthEquation(model, state, state0, dt)
            
            state.widthEq = 1/dt*(state.delta - state0.delta) - state.v;
            
        end
        
        
        function state = updateMassAccumTerm(model, state, state0, dt)

            op = model.operators;
            
            c      = state.c;
            c0     = state0.c;
            delta  = state.delta;
            delta0 = state0.delta;
            
            delta  = op.mapToSei*delta;
            delta0 = op.mapToSei*delta0;
            
            state.massAccum = 1/dt*op.vols.*delta.*(delta.*c - delta0.*c0);
            
        end
            
        function state = updateMassConservation(model, state)
           
            op = model.operators;
            
            flux       = state.flux;
            massSource = state.massSource;
            massAccum  = state.massAccum;
            
            state.massCons = massAccum + op.div(flux) - massSource;

        end
        
        function state = updateSEIgrowthVelocity(model, state)
            
            Mw = model.molecularWeight;
            rho = model.density;
            
            R = state.R;
            
            % Note that R is (for the moment) always negative so that we only have sei layer growth.
            state.v = -0.5*R*(Mw/rho);
            
        end
        
        function state = updateFlux(model, state)
            
            op       = model.operators;
            xi       = op.cellCentroids;
            mapToSei = op.mapToSei;
            faceAver = op.faceAverage;
            
            c     = state.c;
            v     = state.v;
            delta = state.delta;
            
            coef = mapToSei*(delta.*v);
            
            state.flux = op.diffFlux(c) + faceAver*(coef.*(1 - xi).*c);
            
        end

        function state = assembleInterfaceBoundaryEquation(model, state)
        % Boundary equation at the interface side, which relates flux with source term.
            
            op = model.operators;

            c     = state.c;
            cInt  = state.cInterface;
            R     = state.R;
            v     = state.v;
            delta = state.delta;
            
            % Here, xi = 0
            eq = - op.TIntBc.*(op.mapToIntBc*c - cInt) + delta.*v.*cInt - delta.*R;
            
            state.interfaceBoundaryEq = eq;
            
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
    

