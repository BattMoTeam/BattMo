classdef FullSolidDiffusionModel < SolidDiffusionModel

    properties

        %% Input parameters

        % Standard Parameters


        % Function to update diffusion coefficient value, given as a struct with fields
        % - type         : element in {'function', 'constant'}. If 'constant' is chosen the value of referenceDiffusionCoefficient defined in parent class is used
        % - functionname : matlab function name (should be available in path)
        % - argumentlist : should be  ["c", "cmin", "cmax"]
        diffusionCoefficient

        saturationConcentration % the saturation concentration of the guest molecule in the host material (symbol: cmax)
        guestStoichiometry100   % the ratio of the concentration of the guest molecule to the saturation concentration
                                % of the guest molecule in a phase at a cell voltage that is defined as 100% SOC(symbol: theta100)
        guestStoichiometry0     % the ratio of the concentration of the guest molecule to the saturation concentration
                                % of the guest molecule in a phase at a cell voltage that is defined as 0% SOC (symbol: theta0)

        % Advanced parameters

        np             % Number of particles
        N              % Discretization parameters in spherical direction


        %% Computed parameters at model setup

        useDFunc
        computeDFunc % used when useDFunc is true. Function handler to compute D as function of cElectrode, see method updateDiffusionCoefficient


    end

    methods

        function model = FullSolidDiffusionModel(inputparams)

            model = model@SolidDiffusionModel(inputparams);

            fdnames = {'volumeFraction'         , ...
                       'diffusionCoefficient'   , ...
                       'saturationConcentration', ...
                       'guestStoichiometry100'  , ...
                       'guestStoichiometry0'    , ...
                       'np'                     , ...
                       'N'};

            model = dispatchParams(model, inputparams, fdnames);
            model.operators = model.setupOperators();

            D = inputparams.diffusionCoefficient;
            if ~isempty(D)
                switch D.type
                  case 'constant'
                    model.useDFunc = false;
                  case 'function'
                    model.useDFunc = true;
                    model.computeDFunc = str2func(D.functionname);
                  otherwise
                    errror('type of D not recognized.')
                end
            else
                model.useDFunc = false;
            end

        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@SolidDiffusionModel(model);

            varnames = {};
            % concentration in the particle
            varnames{end + 1} = 'c';
            % Average concentration in the particle (not used in assembly)
            varnames{end + 1} = 'cAverage';
            % flux term
            varnames{end + 1} = 'flux';
            % Mass conservation equation
            varnames{end + 1} = 'solidDiffusionEq';

            model = model.registerVarNames(varnames);

            fn = @FullSolidDiffusionModel.updateDiffusionCoefficient;
            if model.useDFunc
                inputnames = {'c'};
            else
                inputnames = {'T'};
            end
            model = model.registerPropFunction({'D', fn, inputnames});

            fn = @FullSolidDiffusionModel.updateFlux;
            inputnames = {'c', 'D'};
            model = model.registerPropFunction({'flux', fn, inputnames});

            fn = @FullSolidDiffusionModel.updateMassConservation;
            inputnames = {'massAccum', 'flux', 'massSource'};
            model = model.registerPropFunction({'massCons', fn, inputnames});

            fn = @FullSolidDiffusionModel.updateMassSource;
            model = model.registerPropFunction({'massSource', fn, {'Rvol'}});

            fn = @FullSolidDiffusionModel.updateMassAccum;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'massAccum', fn, {'c'}});

            fn = @FullSolidDiffusionModel.assembleSolidDiffusionEquation;
            model = model.registerPropFunction({'solidDiffusionEq', fn, {'c', 'cSurface', 'massSource', 'D'}});

            fn = @FullSolidDiffusionModel.updateAverageConcentration;
            model = model.registerPropFunction({'cAverage', fn, {'c'}});

            % we remove this declaration as it is not used in assembly (otherwise it may be computed but not used)
            model = model.setAsExtraVarName('cAverage');

        end

        function operators = setupOperators(model)

            np = model.np;
            N  = model.N;
            rp = model.particleRadius;

            celltbl.cells = (1 : np)';
            celltbl = IndexArray(celltbl);

            % Solid particle cells
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

            radii = G.cells.centroids;

            G.faces.centroids = r;
            G.faces.areas     = 4*pi*r.^2;
            G.faces.normals   = G.faces.areas;

            rock.perm = ones(N, 1);
            rock.poro = ones(N, 1);

            tbls = setupTables(G);
            cellfacetbl = tbls.cellfacetbl;

            hT = computeTrans(G, rock); % hT is in cellfacetbl

            cells = cellfacetbl.get('cells');
            faces = cellfacetbl.get('faces');
            sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1; % sgn is in cellfacetbl

            % We change name of cellfacetbl
            ScellSfacetbl = cellfacetbl;
            ScellSfacetbl = replacefield(ScellSfacetbl, {{'cells', 'Scells'}, {'faces', 'Sfaces'}});

            % Here, we use that we know *apriori* the indexing in G.cells.faces (the last index corresponds to outermost cell-face)
            Tbc = hT(end); % half-transmissibility for of the boundary face
            Tbc = repmat(Tbc, np, 1);

            Sfacetbl.Sfaces = (2 : N)'; % index of the internal faces (correspond to image of C')
            Sfacetbl = IndexArray(Sfacetbl);
            cellSfacetbl = crossIndexArray(celltbl, Sfacetbl, {}, 'optpureproduct', true);

            % we consider only the internal faces
            allScellSfacetbl = ScellSfacetbl;
            ScellSfacetbl = crossIndexArray(allScellSfacetbl, Sfacetbl, {'Sfaces'});

            cellScellSfacetbl = crossIndexArray(celltbl, ScellSfacetbl, {}, 'optpureproduct', true);

            map = TensorMap();
            map.fromTbl = allScellSfacetbl;
            map.toTbl = cellScellSfacetbl;
            map.mergefds = {'Scells', 'Sfaces'};
            map = map.setup();

            hT = map.eval(hT);
            sgn = map.eval(sgn);

            %% setup of divergence operator

            prod = TensorProd();
            prod.tbl1 = cellScellSfacetbl;
            prod.tbl2 = cellSfacetbl;
            prod.tbl3 = cellScelltbl;
            prod.mergefds = {'cells'};
            prod.reducefds = {'Sfaces'};
            prod = prod.setup();

            divMat = SparseTensor();
            divMat = divMat.setFromTensorProd(sgn, prod);
            divMat = divMat.getMatrix();

            div = @(u) (divMat*u);

            gradMat = -divMat';

            prod = TensorProd();
            prod.tbl1 = cellScellSfacetbl;
            prod.tbl2 = cellScelltbl;
            prod.tbl3 = cellSfacetbl;
            prod.mergefds = {'cells'};
            prod.reducefds = {'Scells'};
            prod = prod.setup();

            invHtMat = SparseTensor();
            invHtMat = invHtMat.setFromTensorProd(1./hT, prod);
            invHtMat = invHtMat.getMatrix();

            flux = @(D, c) -(1./(invHtMat*(1./D))).*(gradMat*c);

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

            %% map from cell (celltbl) to cell-particle (cellScelltbl)
            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'cells'};
            map = map.setup();

            mapToParticle = SparseTensor();
            mapToParticle = mapToParticle.setFromTensorMap(map);
            mapToParticle = mapToParticle.getMatrix();

            vols = G.cells.volumes;

            map = TensorMap();
            map.fromTbl = Scelltbl;
            map.toTbl = cellScelltbl;
            map.mergefds = {'Scells'};
            map = map.setup();

            vols = map.eval(vols);

            operators = struct('div'          , div          , ...
                               'flux'         , flux         , ...
                               'mapFromBc'    , mapFromBc    , ...
                               'mapToParticle', mapToParticle, ...
                               'mapToBc'      , mapToBc      , ...
                               'Tbc'          , Tbc          , ...
                               'vols'         , vols         , ...
                               'radii'        , radii);

        end

        function jsonstruct = exportParams(model)

            jsonstruct = exportParams@SolidDiffusionModel(model);

            fdnames = {'diffusionCoefficient'   , ...
                       'saturationConcentration', ...
                       'guestStoichiometry100'  , ...
                       'guestStoichiometry0'    , ...
                       'useDFunc'               , ...
                       'computeDFunc'};

            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                jsonstruct.(fdname) = model.(fdname);
            end

        end


        function state = updateMassSource(model, state)

            op  = model.operators;
            rp  = model.particleRadius;
            vf  = model.volumeFraction;

            Rvol = state.Rvol;

            Rvol = op.mapFromBc*Rvol;

            state.massSource = - Rvol.*((4*pi*rp^3)./(3*vf));

        end

        function state = updateDiffusionCoefficient(model, state)

            if model.useDFunc

                computeD = model.computeDFunc;
                cmax     = model.saturationConcentration;
                theta0   = model.guestStoichiometry0;
                theta100 = model.guestStoichiometry100;

                c = state.c;

                cmin = theta0*cmax;
                cmax = theta100*cmax;

                soc = (c - cmin)./(cmax - cmin);

                D = computeD(soc);

                state.D = D;

            else

                state = updateDiffusionCoefficient@SolidDiffusionModel(model, state);

            end

        end

        function state = updateMassAccum(model, state, state0, dt)

            op = model.operators;

            c = state.c;
            c0 = state0.c;

            state.massAccum = 1/dt*op.vols.*(c - c0);

        end

        function state = updateMassConservation(model, state)

            op = model.operators;

            flux       = state.flux;
            massSource = state.massSource;
            massAccum  = state.massAccum;

            state.massCons = massAccum + op.div(flux) - massSource;

        end

        function state = updateFlux(model, state)

            useDFunc = model.useDFunc;
            op = model.operators;

            c = state.c;
            D = state.D;

            if useDFunc
                state.flux = op.flux(D, c);
            else
                D = op.mapToParticle*D;
                state.flux = op.flux(D, c);
            end


        end

        function state = assembleSolidDiffusionEquation(model, state)

        %% TODO : change name of this function

            op       = model.operators;
            useDFunc = model.useDFunc;

            c     = state.c;
            D     = state.D;
            cSurf = state.cSurface;
            src   = state.massSource;

            if ~useDFunc
                % TODO : make this implementation better
                % Here, we first dispatch D on all the particle cells and, then, retain only the value at the boundary.
                D = op.mapToParticle*D;
            end

            D = op.mapToBc*D;

            eq = D.*op.Tbc.*(op.mapToBc*c - cSurf) + op.mapToBc*src;

            state.solidDiffusionEq = eq;

        end

        function state = updateAverageConcentration(model, state)

            op = model.operators;
            vols = op.vols;
            map = op.mapToParticle;

            c = state.c;

            m    = map'*(c.*vols); % total amount [mol] in the cell particles
            vols = map'*(vols);    % volume

            state.cAverage = m./vols;

        end

        function c = getParticleConcentrations(model, state)
        % Reshape the particle concentration distribution as an array
            np = model.np;
            N  = model.N;

            c = reshape(state.c, np, N)';
            
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
