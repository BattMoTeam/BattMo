classdef MaxwellStefanDiffusion < BaseModel
% Maxwell stefan diffusion model as provided by
% @article{curtiss1999multicomponent,
%   title={Multicomponent diffusion},
%   author={Curtiss, CF and Bird, R Byron},
%   journal={Industrial \& Engineering Chemistry Research},
%   volume={38},
%   number={7},
%   pages={2515--2522},
%   year={1999},
%   publisher={ACS Publications}
% }
%
% We consider (for now) isothermal case and no external forces
    
    properties

        % List of pair of diffusion coefficients. The format is a struct with fields
        % - values 
        % - unit
        % The field values is a list of cells. Each cell in the list is itself a cell, where
        % - cell 1 : name of the first component (string)
        % - cell 2 : name of the second component (string)
        % - cell 3 : value of the binary diffusion coefficient between the two components
        % This structure is processed in setupDiffusionMatrix to obtain the diffusionMatrix 
        diffusionCoefficients 
        % Cell array with the component names
        componentNames
        % Struct which maps a component name (provided as a struct field name) to component index (value returned by the
        % struct)
        compInds
        % Number of components
        numberOfComponents
        % Molar weights (array with N values, where N is the number of components)
        molarWeights
        % Permeability
        permeability
        % Viscosity
        viscosity
        
        % Boundary sub-model
        Boundary

        % Helpers
        
        % Diffusion matrix coefficients of dimension NxN (where N is the number of components)
        diffusionMatrix
        
        boundaryFaces % Boundary faces that are included with their own degrees of freedom. For the boundary faces that are
                      % not included we have no flux boundary conditions

        mappings % Mappings given by the followings fields
                 % - mapToBc   : injection mapping from domain cell value (dimension model.grid.cells.num) to boundary face values
                 %             (dimension numel(boundaryFaces))
                 % - mapFromBc : injection mapping from boundary face values
                 %             (dimension numel(boundaryFaces)) to domain cell value (dimension model.grid.cells.num)
    end
    
    methods

        function model = MaxwellStefanDiffusion(inputparams)

            model = model@BaseModel();
            
            fdnames = {'G'              , ...
                       'diffusionCoefficients', ...
                       'componentNames' , ...
                       'permeability'   , ...
                       'viscosity'      , ...
                       'molarWeights'};
            
            model = dispatchParams(model, inputparams, fdnames);

            model.numberOfComponents = numel(model.componentNames);

            model.subModelNameList = {'Boundary'};
            
            ncomp = model.numberOfComponents;
            
            for icomp = 1 : ncomp

                name = model.componentNames{icomp};
                compInds.(name) = icomp;
                
            end

            model.compInds = compInds;
            
            model = setupDiffusionMatrix(model);
            
        end


        function model = setupDiffusionMatrix(model)

            dcs      = model.diffusionCoefficients.values;
            ncomp    = model.numberOfComponents;
            compinds = model.compInds;
            u        = eval(model.diffusionCoefficients.unit);

            D = nan(ncomp, ncomp);
            for idc = 1 : numel(dcs)
                dc = dcs{idc};
                i1 = compinds.(dc{1});
                i2 = compinds.(dc{2});
                v = dc{3}*u;
                D(i1, i2) = v;
                D(i2, i1) = v;
            end

            model.diffusionMatrix = D;
            
        end
        
        function model = setupMappings(model)

        %  boundaryFaces gives the indexing of the boundary faces

            bcfacetbl.faces = model.boundaryFaces;
            bcfacetbl = IndexArray(bcfacetbl);
            bcfacetbl = bcfacetbl.addLocInd('ind'); % add local indexing (used below)
           
            tbls = setupTables(model.grid);

            bccellfacetbl = crossIndexArray(bcfacetbl, tbls.cellfacetbl, {'faces'});
            bccellfacetbl = sortIndexArray(bccellfacetbl, {'ind'}, 'keepAllFields', true);

            celltbl = tbls.celltbl;

            map = TensorMap();
            map.fromTbl  = bccellfacetbl;
            map.toTbl    = celltbl;
            map.mergefds = {'cells'};
            map = map.setup();

            mapFromBc = map.getMatrix();

            map = TensorMap();
            map.fromTbl  = celltbl;
            map.toTbl    = bccellfacetbl;
            map.mergefds = {'cells'};
            map = map.setup();

            mapToBc = map.getMatrix();

            model.mappings = struct('mapFromBc', mapFromBc, ...
                                    'mapToBc'  , mapToBc);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            ncomp = model.numberOfComponents;

            % Define first the variables that are common to the main domain and the boundary
            varnames = {};
            % Chemical potential
            varnames{end + 1} = VarName({}, 'chemicalPotentials', ncomp);
            % Diffusional driving forces
            varnames{end + 1} = VarName({}, 'diffusionForces', ncomp);            
            % Diffusion fluxes
            varnames{end + 1} = VarName({}, 'diffusionFluxes', ncomp);
            % Convection fluxes
            varnames{end + 1} = VarName({}, 'convectionFluxes', ncomp);
            % Mass fluxes (sum over all fluxes)
            varnames{end + 1} = VarName({}, 'massFluxes', ncomp);
            % Maxwell Stefan equation. They relate the driving forces with the diffusion fluxes
            varnames{end + 1} = VarName({}, 'maxwellStefanEquations', ncomp);
            % Mass average velocity vector
            varnames{end + 1} = 'velocity';
            % Mass fractions
            varnames{end + 1} = VarName({}, 'massFractions', ncomp);
            % Mass fraction equation. The mass fractions sum to one.
            varnames{end + 1} = 'massFractionEquation';
            % Pressure / Pa
            varnames{end + 1} = 'pressure';
            % Temperature / K
            varnames{end + 1} = 'temperature';
            % Total density / kg/m^3
            varnames{end + 1} = 'density';            
            % Concentration / mol/m^3
            varnames{end + 1} = 'concentration';
            % Total molecular weight / kg/mol
            varnames{end + 1} = 'molarWeight';
            
            model = model.registerVarNames(varnames);
            
            for ivar = 1 : numel(varnames)
                varname = varnames{ivar};
                if ~isa(varname, 'VarName')
                    varname = VarName({'Boundary'}, varname);
                else
                    varname.namespace = {'Boundary'};
                end
                varnames{ivar} = varname;
            end
            
            model = model.registerVarNames(varnames);
            
            % Register the variable specific to main domain (mass conservation equations)
            varnames = {};
            % Mass accumulation terms
            varnames{end + 1} = VarName({}, 'massAccums', ncomp);
            % Mass source terms
            varnames{end + 1} = VarName({}, 'sources', ncomp);
            % Boundary mass source terms
            varnames{end + 1} = VarName({}, 'boundarySources', ncomp);
            % Mass conservation equations
            varnames{end + 1} = VarName({}, 'massConses', ncomp);
            
            model = model.registerVarNames(varnames);
            
            model = model.setAsStaticVarName('temperature');
            model = model.setAsStaticVarName({'Boundary', 'temperature'});

            fn = @MaxwellStefanGasDiffusion.updateMolarWeight;
            inputvarnames = {VarName({}, 'massFractions', ncomp)};
            model = model.registerPropFunction({'molarWeight', fn, inputvarnames});

            fn = @MaxwellStefanGasDiffusion.updateMolarWeight;
            inputvarnames = {VarName({'Boundary'}, 'massFractions', ncomp)};
            model = model.registerPropFunction({{'Boundary', 'molarWeight'}, fn, inputvarnames});

            fn = @MaxwellStefanGasDiffusion.updateVelocity;
            inputvarnames = {VarName({}, 'pressure')};
            model = model.registerPropFunction({'velocity', fn, inputvarnames});

            fn = @MaxwellStefanGasDiffusion.updateVelocity;
            inputvarnames = {VarName({'Boundary'}, 'pressure')};
            model = model.registerPropFunction({{'Boundary', 'velocity'}, fn, inputvarnames});
            
            fn = @MaxwellStefanDiffusion.updateConcentration;
            model = model.registerPropFunction({'concentration', fn, {'density', 'molarWeight'}});
            model = model.registerPropFunction({{'Boundary', 'concentration'}, fn, {{'Boundary', 'density'}, {'Boundary', 'molarWeight'}}});
            
            fn = @MaxwellStefanDiffusion.updateMassFractionEquation;
            inputvarnames = {VarName({}, 'massFractions', ncomp)};
            model = model.registerPropFunction({'massFractionEquation', fn, inputvarnames});

            fn = @MaxwellStefanDiffusion.updateMassFractionEquation;
            inputvarnames = {VarName({'Boundary'}, 'massFractions', ncomp)};
            outputvarname = VarName({'Boundary'}, 'massFractionEquation');
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @MaxwellStefanDiffusion.updateMaxwellStefanEquations;
            inputvarnames = {VarName({}, 'diffusionForces', ncomp), ...
                             VarName({}, 'diffusionFluxes', ncomp), ...
                             VarName({}, 'molarWeight')           , ...                             
                             VarName({}, 'massFractions', ncomp)};
            outputvarname = VarName({}, 'maxwellStefanEquations', ncomp);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @MaxwellStefanDiffusion.updateMaxwellStefanEquations;
            inputvarnames = {VarName({'Boundary'}, 'diffusionForces', ncomp), ...
                             VarName({'Boundary'}, 'diffusionFluxes', ncomp), ...
                             VarName({'Boundary'}, 'molarWeight')           , ...                             
                             VarName({'Boundary'}, 'massFractions', ncomp)};
            outputvarname = VarName({'Boundary'}, 'maxwellStefanEquations', ncomp);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            
            for icomp = 1 : ncomp
                                
                fn = @MaxwellStefanDiffusion.updateConvectionFluxes;
                outputvarname = VarName({}, 'convectionFluxes', ncomp, icomp);
                inputvarnames  = {'density'                                  , ...
                                  'velocity', ...
                                  VarName({}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanDiffusion.updateConvectionFluxes;
                outputvarname = VarName({'Boundary'}, 'convectionFluxes', ncomp, icomp);
                inputvarnames  = {VarName({'Boundary'},'density') , ...
                                  VarName({'Boundary'},'velocity'), ...
                                  VarName({'Boundary'}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanDiffusion.updateMassFluxes;
                outputvarname = VarName({}, 'massFluxes', ncomp, icomp);
                inputvarnames  = {VarName({}, 'convectionFluxes', ncomp, icomp), ...
                                  VarName({}, 'diffusionFluxes', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});                

                fn = @MaxwellStefanDiffusion.updateMassFluxes;
                outputvarname = VarName({'Boundary'}, 'massFluxes', ncomp, icomp);
                inputvarnames  = {VarName({'Boundary'}, 'convectionFluxes', ncomp, icomp), ...
                                  VarName({'Boundary'}, 'diffusionFluxes', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});                
                
                fn = @MaxwellStefanDiffusion.updateDiffusionForces;
                outputvarname = VarName({}, 'diffusionForces', ncomp, icomp);
                inputvarnames  = {'temperature'                                  , ...
                                  'pressure'                                     , ...
                                  'molarWeight'                                  , ...
                                  'concentration'                                , ...
                                  VarName({}, 'chemicalPotentials', ncomp, icomp), ...
                                  VarName({}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanDiffusion.updateBoundaryDiffusionForces;
                outputvarname = VarName({'Boundary'}, 'diffusionForces', ncomp, icomp);
                inputvarnames  = {'temperature'                                            , ...
                                  'pressure'                                               , ...
                                  'density'                                                , ...
                                  'concentration'                                          , ...
                                  VarName({}, 'chemicalPotentials', ncomp, icomp)          , ...
                                  VarName({}, 'massFractions', ncomp, icomp)               , ...
                                  VarName({'Boundary'}, 'temperature')                     , ...
                                  VarName({'Boundary'}, 'pressure')                        , ...
                                  VarName({'Boundary'}, 'density')                         , ...
                                  VarName({'Boundary'}, 'concentration')                   , ...
                                  VarName({'Boundary'}, 'chemicalPotentials', ncomp, icomp), ...
                                  VarName({'Boundary'}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @MaxwellStefanDiffusion.updateMassAccums;
                fn = {fn, @(propfunc) PropFunction.accumFuncCallSetupFn(propfunc)};
                outputvarname = VarName({}, 'massAccums', ncomp, icomp);
                inputvarnames  = { 'density', ...
                                   VarName({}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanDiffusion.updateBoundarySources;
                outputvarname = VarName({}, 'boundarySources', ncomp, icomp);
                inputvarnames  = {VarName({'Boundary'}, 'diffusionFluxes', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanDiffusion.updateMassConses;
                outputvarname = VarName({}  , 'massConses'     , ncomp, icomp);
                inputvarnames  = {VarName({}, 'massAccums'     , ncomp, icomp), ...
                                  VarName({}, 'massFluxes'     , ncomp, icomp), ...
                                  VarName({}, 'sources'        , ncomp, icomp), ...
                                  VarName({}, 'boundarySources', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});                
                
            end
            
        end

        function state = updateMaxwellStefanEquations(model, state)

            D   = model.diffusionMatrix;
            mws = model.molarWeights;

            function eqs = setupEquations(mw, mfs, rho, jds, dfs)
            % mw  = molarWeight;
            % mfs = massFractions;
            % rho = density;
            % jds = diffusionFluxes;
            % dfs = diffusionForces;
                
                coef = mw^2./rho;
                
                for icomp = 1 : ncomp

                    eqs{icomp} = dfs{icomp};

                    for jcomp = 1 : ncomp

                        if jcomp ~= icomp

                            eqs{icomp} = eqs{icomp} -  ...
                                coef./D(icomp, jcomp)*( ...
                                mfs{jcomp}.*(mws(jcomp)/mws(icomp)).*jds{icomp}  ...
                                - mfs{icomp}.*(mws(icomp)/mws(jcomp)).*jds{jcomp});
                            
                        end

                    end

                end

                % We replace the first of thes equations, which otherwise are not independent due to symmetry of D(i,j),
                % with the null average diffusion flux equation

                eqs{1} = jds{1};

                for icomp = 2 : ncomp
                    eqs{1} = eqs{1} + jds{icomp};
                end

            end
            
            mw  = state.molarWeight;
            mfs = state.massFractions;
            rho = state.density;
            jds = state.diffusionFluxes;
            dfs = state.diffusionForces;
            
            state.maxwellStefanEquations = setupEquations(mw, mfs, rho, jds, dfs);

            mw  = state.Boundary.molarWeight;
            mfs = state.Boundary.massFractions;
            rho = state.Boundary.density;
            jds = state.Boundary.diffusionFluxes;
            dfs = state.Boundary.diffusionForces;

            state.Boundary.maxwellStefanEquations = setupEquations(mw, mfs, rho, jds, dfs);
            
        end
        
        function state = updateMassFluxes(model, state)

            for icomp = 1 : ncomp

                state.massFluxes{icomp}          = state.convectionFluxes{icomp} + state.diffusionFluxes{icomp};
                state.Boundary.massFluxes{icomp} = state.Boundary.convectionFluxes{icomp} + state.Boundary.diffusionFluxes{icomp};
                
            end
            
        end

        
        function state = updateConvectionFluxes(model, state)

            ncomp = model.numberOfComponents();
            
            for icomp = 1 : ncomp

                state.convectionFluxes{icomp}          = state.massFractions{icomp}.*state.density.*state.velocity;
                state.Boundary.convectionFluxes{icomp} = state.Boundary.massFractions{icomp}.*state.Boundary.density.*state.Boundary.velocity;
                
            end
            
        end
            
        function state = updateVelocity(model, state)

            K       = model.permeability;
            mu      = model.viscosity;
            bcfaces = model.boundaryFaces;
            
            p    = state.pressure
            bd_p = state.Boundary.pressure;

            mob = K/mu; % mobility

            state.velocity          = assembleFlux(model, p, mob);
            state.Boundary.velocity = assembleFlux(model, p, bd_p, mob, bcfaces);
            
        end

        
        function state = updateMassAccums(model, state, state0, dt)

            ncomp = model.numberOfComponents;

            for icomp = 1 : ncomp

                state.massAccums{icomp} = 1/dt*(state.density.*state.massFractions{icomp} - ...
                                                state0.density.*state0.massFractions{icomp});
            end
            
        end
        
        function state = updateMassFractionEquation(model, state)

            ncomp = model.numberOfComponents;
            
            state.massFractionEquation          = 1;
            state.Boundary.massFractionEquation = 1;

            for icomp = 1 : ncomp
                state.massFractionEquation          = state.massFractionEquation          - state.massFractions{icomp};
                state.Boundary.massFractionEquation = state.Boundary.massFractionEquation - state.Boundary.massFractions{icomp};
            end
            
        end

        function state = updateMolarWeight(model, state)

            state.molarWeight          = model.computeMolarWeight(state.massFractions);
            state.Boundary.molarWeight = model.computeMolarWeight(state.Boundary.massFractions);
            
        end

        function molarWeight = computeMolarWeight(model, massFractions)

            ncomp = model.numberOfComponents;
            mws   = model.molarWeights;
            
            mw = 0*mfs{1}; % (useful?) trick to make sure we have AD instantiation when needed.

            for icomp = 1 : ncomp

                mw = mw + massFractions{icomp}/mws(icomp);
                
            end
            
            molarWeight = 1./mw;
            
        end

        function state = updateConcentration(model, state)

            state.concentration          = state.density./state.molarWeight;
            state.Boundary.concentration = state.Boundary.density./state.Boundary.molarWeight;
            
        end

        function state = updateDiffusionFluxes(model, state)

            
            
        end
        
        function state = updateDiffusionForces(model, state)
        % Equation 1.9 in main reference (Curtiss and Bird, 1999)

            R = PhysicalConstants.R;
            
            mws = model.molarWeights;
            
            mw  = state.molarWeight;
            T   = state.T;
            c   = state.concentration;
            mfs = state.massFractions;
            p   = state.pressure;
            mus = state.chemicalPotentials;
            
            for icomp = 1 : ncomp

                potential        = 1./T.*mus{icomp}/mws(icomp);
                fluxCoeffficient = mfs{icomp}.*mw/R;

                % Note minus sign, because in there is already a minus sign in the assembleFlux function
                dcomp = - assembleFlux(model, potential, fluxCoefficient);

                potential        = p;
                fluxCoeffficient = mfs{icomp}./(R*c.*T);
                
                dcomp = dcomp + assembleFlux(model, potential, fluxCoefficient);

                state.diffusionForces{icomp} = dcomp;
                
            end

        end

        function state = updateBoundaryDiffusionForces(model, state)

        % We set the diffustionForces at the boundary as an equation. Otherwise, the implementation is consistent with
        % updateDiffusionForces in the interior of the domain
            
            R = PhysicalConstants.R;

            ncomp   = model.numberOfComponents;
            bcfaces = model.boundaryFaces;
            mws     = model.molarWeights;

            nbf = numel(bcfaces);

            % short names for variables (interior)
            mw  = state.molarWeight;
            T   = state.T;
            c   = state.concentration;
            mfs = state.massFractions;
            p   = state.pressure;
            mus = state.chemicalPotentials;
            
            % short names for variables (boundary)
            bd_T   = state.Boundary.T;
            bd_mus = state.Boundary.chemicalPotentials;            
            bd_mfs = state.Boundary.massFractions;

            zeroAD = model.AutoDiffBackend.convertToAD(zeros(nbf, 1), p);
            
            for icomp = 1 : ncomp
                
                dforces{icomp} = zeroAD;

                potential        = 1./T.*mus{icomp}/mws(icomp);
                bd_potential     = 1./bd_T.*bd_mus{icomp}/mws(icomp);
                fluxCoeffficient = mfs{icomp}.*mw/R;

                % Note minus sign, because in there is already a minus sign in the assembleBoundaryFlux function
                dforces{icomp} = dforces{icomp} - assembleBoundaryFlux(potential, bd_potential, fluxCoefficient, bcfaces);

                potential        = p;
                bd_potential     = bd_p;
                fluxCoeffficient = mfs{icomp}./(R*c.*T);
                
                dforces{icomp} = dforces{icomp} + assembleBoundaryFlux(potential, bd_potential, fluxCoefficient, bcfaces);

            end
            
            state.diffusionForces = dforces;

        end

        function state = updateBoundarySources(model, state)

            mapFromBc = model.mappings.mapFromBc;
            ncomp     = model.numberOfComponents;
            
            for icomp = 1 : ncomp

                % note the minus sign as the boundary flux are computed outwards (from interior to outer), see function assembleBoundaryFlux.
                state.boundaryForces{icomp} = - mapFromBc*state.Boundary.diffusionFluxes{icomp};

            end
            
        end
        
    end

end

    

