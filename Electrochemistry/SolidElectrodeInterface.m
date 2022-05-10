classdef SolidElectrodeInterface < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();
        
        molecularWeight % SEI molecular weight [kg/mol]
        density         % SEI densisity [kg/m^3]
        D               % SEI diffusion coefficient [m^2/s]
        cExternal       % SEI concentraton at the outer side
        
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
                       'D'};
            model = dispatchParams(model, paramobj, fdnames);
            
            %% FIXME : fix the operators
            % model.operators = model.setupOperators();
            
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
            % surface concentration
            varnames{end + 1} = 'seiwidth';
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
            % Mass conservation equation
            varnames{end + 1} = 'solidDiffusionEq';
            % SEI Width equation
            varnames{end + 1} = 'widthEq';
            
            model = model.registerVarNames(varnames);

            fn = @SolidElectrodeInterface.updateFlux;
            inputnames = {'c', 'v'};
            model = model.registerPropFunction({'flux', fn, inputnames});
            
            fn = @SolidElectrodeInterface.updateSEIgrowthVelocity;
            inputnames = {'R'};
            model = model.registerPropFunction({'v', fn, inputnames});
            
            fn = @SolidElectrodeInterface.updateMassConservation;
            inputnames = {'massAccum', 'flux', 'massSource'};
            model = model.registerPropFunction({'massCons', fn, inputnames});

            fn = @SolidElectrodeInterface.updateMassSource;
            model = model.registerPropFunction({'massSource', fn, {'R', 'c', 'v'}});
            
            fn = @SolidElectrodeInterface.assembleSolidDiffusionEquation;
            model = model.registerPropFunction({'solidDiffusionEq', fn, {'c', 'cSurface', 'massSource'}});
            
            fn = @SolidElectrodeInterface.assembleWidthEquation;
            model = model.registerPropFunction({'widthEq', fn, {'seiwidth', 'v'}});
        end
        
        function operators = setupOperators(model)
            
        % FIXME : fix operators
            
            np = model.np;
            N  = model.N;
            
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
            
            rock.perm = ones(N, 1);
            rock.poro = ones(N, 1);

            op = setupOperatorsTPFA(G, rock);
            C = op.C;
            T = op.T; % in Sfacetbl
            T_all = op.T_all;
            
            TExtBc = T_all(N); % half-transmissibility for of the external boundary face
            TIntBc = T_all(1); % half-transmissibility for of the internal boundary face (close to solid particle)
            
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

            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = cellSfacetbl;
            map.mergefds = {'cells'};
            map = map.setup();

            flux = @(D, delta, c) - map.eval(D./delta).*(Grad*c);

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
            
            operators = struct('div'         , div         , ...
                               'flux'        , flux        , ...
                               'mapFromExtBc', mapFromExtBc, ...
                               'mapToExtBc'  , mapToExtBc  , ...
                               'mapFromIntBc', mapFromIntBc, ...
                               'mapToIntBc'  , mapToIntBc  , ...
                               'TIntBc'      , TIntBc      , ...
                               'TExtBc'      , TExtBc);
            
        end

        function state = updateMassSource(model, state)
            
            op = model.operators;
            rp = model.rp;
            volumetricSurfaceArea = model.volumetricSurfaceArea;
            cExternal = model.cExternal;
            
            R = state.R;
            c = state.c;
            
            c = op.mapToExtBc*c;
            srcExternal = op.TExtBc.*(c - cExternal) + 0.5*v.*(c + cExternal);
            srcExternal = op.mapFromExtBc*srcExternal;
            
            %% FIXME : check expression
            srcInterface = -op.mapFromIntBc*R;
            state.massSource = srcExternal + srcInterface;
            
        end

        
        function state = assembleWidthEquation(model, state, state0, dt)
            
            state.widthEq = 1/dt*(state.seiwidth - state0.seiwidth) - state.v
            
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
        
        function state = updateSEIgrowthVelocity(model, state)
            
            Mw = model.molecularWeight;
            rho = model.density;
            
            R = state.R;
            
            state.v = -0.5*R*(Mw/rho);
            
        end
        
        
        function state = updateFlux(model, state)
            
            op = model.operators;
            D = model.D;
            
            c = state.c;
            v = state.v;
            
            % FIXME : check op.average exists (mapping from cell to internal face)
            % FIXME : map v from sidereaction "space" to mass cons "space"
            state.flux = op.flux(D, c) - op.average(v.*c);
            
        end

        function state = assembleSolidDiffusionEquation(model, state)
            
            op = model.operators;

            c     = state.c;
            cSurf = state.cSurface;
            src   = state.massSource;
            v     = state.v;
            
            % FIXME : check sign before v, use average of concentration to compute convection term (?)
            eq = op.TIntBc.*(op.mapToIntBc*c - cSurf) + v.*cSurf + op.mapToIntBc*src;
            
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
    

