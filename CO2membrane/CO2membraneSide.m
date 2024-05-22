classdef CO2membraneSide < BaseModel

    properties

        Boundary

        couplingTerms
        
        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        viscosity   % viscosity
        temperature % temperature
        diameter    % tube parameter

        control % Structure with fields
                % - pressure
                % - rate
                % - composition
        
        % Advanced parameter
        poiseuilleCoefficient

        %% Helpers
        
        boundaryHelper % control stucture with fields
                       % - transmissibilities
                       % - mapToBc
                       % - mapFromBc
                       % - molFractions (Vector of dimension nGas with the mol fractions. It used for the moment to setup
                       %                 the composition at all boundaries. Fragile!)
        controlHelper % control structure with fields
                      % - pressureMap
                      % - pressureValues
                      % - fluxMap
                      % - fluxValues
    end

    methods
        
        function model = CO2membraneSide(inputparams)
            
            model = model@BaseModel();

            fdnames = {'G'           , ...
                       'couplingTerms', ...
                       'viscosity'   , ...
                       'temperature' , ...
                       'diameter'    , ...
                       'control'     , ...
                       'poiseuilleCoefficient'};

            model = dispatchParams(model, inputparams, fdnames);

            model.constants = PhysicalConstants();

            model.Boundary = CO2membraneSideBoundary([]);
            
            model = CO2membrane.setupGasStructures(model);
            
            model = model.setupHelpers();
            
            if isempty(model.poiseuilleCoefficient)
                % setup function to update it
                error('not yet implemented');
            end

            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            nGas = model.nGas;
            
            varnames = {};
            %  fluxes [mol/s]
            varnames{end + 1} = VarName({}, 'fluxes', nGas);
            %  mol fractions
            varnames{end + 1} = VarName({}, 'molFractions', nGas);
            %  mol fractions constraint
            varnames{end + 1} = 'molFractionConstraint';
            %  pressures [pascal]
            varnames{end + 1} = VarName({}, 'pressures', nGas);
            %  total pressure [pascal]
            varnames{end + 1} = 'pressure';
            %  boundary sources [mol/s]
            varnames{end + 1} = VarName({}, 'bcSources', nGas);

            % mass conservation equations
            varnames{end + 1} = VarName({}, 'massConses', nGas);
            
            % Boundary mol fraction equation setup
            varnames{end + 1} = VarName({'Boundary'}, 'bcMolFractionDefinitions', nGas);
            
            % Boundary flux definition
            varnames{end + 1} = {'Boundary', 'bcFluxDefinition'};
            
            % Control 
            varnames{end + 1} = VarName({'Boundary'}, 'control');
            
            model = model.registerVarNames(varnames);

            for igas = 1 : nGas
                
                fn = @CO2membraneSide.updateFluxes;
                inputvarnames = {'pressure', VarName({}, 'molFractions', nGas, igas)};
                outputvarname = VarName({}, 'fluxes', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @(model, state) CO2membraneSide.updatePressures(model, state);
                fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
                inputvarnames = {'pressure', VarName({}, 'molFractions', nGas, igas)};
                outputvarname = VarName({}, 'pressures', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @CO2membraneSide.updateSources;
                inputvarnames = {VarName({'Boundary'}, 'fluxes', nGas, igas)};
                outputvarname = VarName({}, 'bcSources', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @CO2membraneSide.updateBoundaryMolFractionsDefinitions;
                inputvarnames = {VarName({'Boundary'}, 'molFractions', nGas, igas), ...
                                 VarName({}, 'molFractions', nGas, igas)          , ...
                                 'pressure'                                       , ...
                                 {'Boundary', 'pressure'}};
                outputvarname = VarName({'Boundary'}, 'bcMolFractionDefinitions', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                if model.isRootSimulationModel

                    fn = @CO2membraneSide.updateMassConses;
                    fn = {fn, @(prop) PropFunction.accumFuncCallSetupFn(prop)};
                    inputvarnames = {VarName({}, 'molFractions', nGas, igas)   , ...
                                     VarName({}, 'fluxes', nGas, igas)   , ...
                                     VarName({}, 'bcSources', nGas, igas)};
                    outputvarname = VarName({}, 'massConses', nGas, igas);
                    model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                end

            end

            fn = @CO2membraneSide.updateMolFractionConstraint;
            intputvarnames = {VarName({}, 'fluxes', nGas)};
            outputvarname = 'molFractionConstraint';
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});            
            
            fn = @CO2membraneSide.updateBcFluxDefinition;
            inputvarnames = {{'Boundary', 'flux'}    , ...
                             {'Boundary', 'pressure'}, ...
                             'pressure'};
            outputvarname = {'Boundary', 'bcFluxDefinition'};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @CO2membraneSide.updateControl;
            inputvarnames = {VarName({'Boundary'}, 'flux'), ...
                             VarName({'Boundary'}, 'pressure')};
            outputvarname = VarName({'Boundary'}, 'control');
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            if model.isRootSimulationModel

                model = model.setAsExtraVarName(VarName({}, 'pressures', nGas));
                
            end

        end

        function model = setupHelpers(model, state)

            gasInd = model.gasInd;
            nGas   = model.nGas;
            
            coupterms = model.couplingTerms;
            coupnames = cellfun(@(coupterm) coupterm.name, coupterms, 'uniformoutput', false);


            ind = strcmp('boundary faces', coupnames);
            coupterm = coupterms{ind};
            
            G = model.G;
            bcfaces = coupterm.couplingfaces;

            Tbc = model.G.getBcTrans(bcfaces);

            % setup the given boundary mol fractions
            fdnames = fieldnames(model.control.composition);
            for igas = 1 : nGas
                fdname = fdnames{igas};
                mfs(gasInd.(fdname)) = model.control.composition.(fdname);
            end
            % We normalize to one
            mfs = mfs./sum(mfs);

            bcfacetbl.faces = bcfaces;
            bcfacetbl = IndexArray(bcfacetbl);

            bccelltbl.cells = coupterm.couplingcells;
            bccelltbl = IndexArray(bccelltbl);
            
            celltbl.cells = (1 : model.G.getNumberOfCells())';
            celltbl = IndexArray(celltbl);

            map = TensorMap();
            map.fromTbl = celltbl;
            map.toTbl = bccelltbl;
            map.mergefds = {'cells'};
            map = map.setup();

            mapToBc = SparseTensor();
            mapToBc = mapToBc.setFromTensorMap(map);
            mapToBc = mapToBc.getMatrix();

            map = TensorMap();
            map.fromTbl = bccelltbl;
            map.toTbl = celltbl;
            map.mergefds = {'cells'};
            map = map.setup();

            mapFromBc = SparseTensor();
            mapFromBc = mapFromBc.setFromTensorMap(map);
            mapFromBc = mapFromBc.getMatrix();
            
            boundaryHelper = struct('transmissibilities', Tbc, ...
                                    'mapFromBc', mapFromBc   , ...
                                    'mapToBc', mapToBc       , ...
                                    'molFractions', mfs);
            
            ind = strcmp('control faces', coupnames);
            coupterm = coupterms{ind};            

            ctrlfacetbl.faces = coupterm.couplingfaces;
            ctrlfacetbl = IndexArray(ctrlfacetbl);

            map = TensorMap();
            map.fromTbl  = bcfacetbl;
            map.toTbl    = ctrlfacetbl;
            map.mergefds = {'faces'};
            map = map.setup();

            M = SparseTensor();
            M = M.setFromTensorMap(map);
            M = M.getMatrix();

            controlHelper = struct('pressureMap', M                        , ...
                                   'fluxMap', M                            , ...
                                   'pressureValues', model.control.pressure, ...
                                   'fluxValues', model.control.rate);
            
            model.boundaryHelper = boundaryHelper;
            model.controlHelper  = controlHelper;
            
        end

        function initstate = setupInitialState(model)
        % Here we give values to the primary variables

            bd = 'Boundary';
            
            nGas     = model.nGas;
            ctrlhelp = model.controlHelper;
            bdhelp   = model.boundaryHelper;
            
            nc  = model.G.getNumberOfCells();
            nbc = numel(model.boundaryHelper.transmissibilities);

            mfs      = bdhelp.molFractions;
            pressure = ctrlhelp.pressureValues(1); % we take the first value as a reasonable guess (in case several were given)
            flux     = ctrlhelp.fluxValues(1);     % we take the first value as a reasonable guess
            
            initstate.pressure = pressure*ones(nc, 1);
            
            initstate.(bd).pressure = pressure*ones(nbc, 1);
            initstate.(bd).flux     = flux*ones(nbc, 1);
            
            for igas = 1 : nGas
                initstate.(bd).molFractions{igas} = mfs(igas)*ones(nbc, 1);
                initstate.molFractions{igas} = mfs(igas)*ones(nc, 1);
            end    
            
        end

        function state = updateMassConses(model, state, state0, dt)

            nGas = model.nGas;
            G    = model.G;
            vols = G.getVolumes();
            
            for igas = 1 : nGas

                
                src = state.bcSources{igas};
                q   = state.fluxes{igas};
                mf  = state.molFractions{igas};
                mf0 = state0.molFractions{igas};
                
                eqs{igas} = vols.*(mf - mf0)/dt + G.getDiv(q) - src;
                
            end

            state.massConses = eqs;

        end

        function state = updateMolFractionConstraint(model, state)

            nGas = model.nGas;

            eq = 1;
            
            for igas = 1 : nGas
                eq = eq - state.molFractions{igas};
            end

            state.molFractionConstraint = eq;
            
        end
        
        function state = updateBcFluxDefinition(model, state)

            nGas    = model.nGas;
            pcoef   = model.poiseuilleCoefficient;
            Tbc     = model.boundaryHelper.transmissibilities;
            mapToBc = model.boundaryHelper.mapToBc;
            
            bd = 'Boundary';
            
            qBc = state.(bd).flux;
            pBc = state.(bd).pressure;
            
            p = state.pressure;
            p = mapToBc*p;

            % By convention, qBc denotes outward flux
            state.(bd).bcFluxDefinition = qBc - pcoef*Tbc.*(p - pBc);
            
        end

        function state = updatePressures(model, state)

            
            nGas = model.nGas;
            
            p = state.pressure;
            mfs = state.volumeFractions;
            
            for igas = 1 : nGas

                pressures{igas} = mfs{igas}.*p;
                
            end

            state.pressures = pressures;

        end
        
        function state = updateSources(model, state)

            bd = 'Boundary';
            
            nGas      = model.nGas;
            mapFromBc = model.boundaryHelper.mapFromBc;
            
            for igas = 1 : nGas

                qbc = state.(bd).fluxes{igas};
                % By convention, qBc are outward boundary fluxes
                ms{igas} = - mapFromBc*qbc;
                
            end

            state.bcSources = ms;

        end

        function state = updateBoundaryMolFractionsDefinitions(model, state)

            bd        = 'Boundary';
            nGas      = model.nGas;
            givenMfBc = model.boundaryHelper.molFractions;
            mapToBc   = model.boundaryHelper.mapToBc;
            
            pBc = state.(bd).pressure;
            p = state.pressure;            
            p = mapToBc*p;

            ind = find(pBc - p >= 0);
            
            for igas = 1 : nGas

                mf = state.molFractions{igas};
                mf = mapToBc*mf;
                
                mfBc = state.(bd).molFractions{igas};

                bcdefs{igas}      = mfBc - mf;
                bcdefs{igas}(ind) = mfBc(ind) - givenMfBc(igas)*ones(numel(ind), 1);
 
            end

            state.(bd).bcMolFractionDefinitions = bcdefs;
                
        end
        
        function state = updateFluxes(model, state)
            
            nGas  = model.nGas;
            pcoef = model.poiseuilleCoefficient;

            p = state.pressure;
            
            for igas = 1 : nGas

                mf = state.molFractions{igas};

                q = assembleHomogeneousFlux(model, p, pcoef);

                qs{igas} = assembleUpwindFlux(model, q, mf);
                
            end
            
            state.fluxes = qs;
        end

        function state = updateControl(model, state)

            bd = 'Boundary';

            ctrlhelp = model.controlHelper;
            
            pBc = state.(bd).pressure;
            qBc = state.(bd).flux;

            eqs{1} = ctrlhelp.pressureMap*pBc - ctrlhelp.pressureValues;
            eqs{2} = ctrlhelp.fluxMap*qBc - ctrlhelp.fluxValues;

            state.(bd).control = vertcat(eqs{:});
            
        end
        
        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end
                
    end

end
    
