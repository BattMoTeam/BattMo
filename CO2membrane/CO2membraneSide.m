classdef CO2membraneSide < BaseModel

    properties

        Boundary

        couplingTerm
        
        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        viscosity   % viscosity
        temperature % temperature
        diameter    % tube parameter
        
        % Advanced parameter
        poiseuilleCoefficient

        %% helpers
        Tbc
        boundarySetup
    end

    methods
        
        function model = CO2membraneSide(inputparams)
            
            model = model@BaseModel();

            fdnames = {'couplingTerm', ...
                       'viscosity'   , ...
                       'temperature' , ...
                       'diameter'    , ...
                       'poiseuilleCoefficient'};

            model = dispatchParams(model, inputparams, fdnames);

            model.constants = PhysicalConstants();

            model.Boundary = CO2membraneSideBoundary([]);
            
            model.gasInd.CO2 = 1;
            model.gasInd.O2  = 2;
            model.gasInd.N2  = 3;
            model.gasInd.Ar  = 4;
            model.nGas = 4;


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

                fn = @(model, state) CO2membraneSide.updateMolFractions(model, state);
                fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
                inputvarnames = {'pressure', VarName({}, 'pressures', nGas, igas)};
                outputvarname = VarName({}, 'molFractions', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @CO2membraneSide.updateSources;
                inputvarnames = {VarName({'Boundary'}, 'fluxes', nGas, igas)};
                outputvarname = VarName({}, 'bcSources', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @CO2membraneSide.updateBoundaryMolFractionsDefinitions;
                inputvarnames = {VarName({}, 'molFractions', nGas, igas), ...
                                 'pressure', ...
                                 {'Boundary', 'pressure'}};
                outputvarname = VarName({'Boundary'}, 'bcMolFractionDefinitions', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            end

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
                
            fn = @(model, state) CO2membraneSide.updatePressure(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            inputvarnames = {VarName({}, 'pressures', nGas)};
            outputvarname = 'pressure';
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
        end

        function state = updateBcFluxDefinition(model, state)

            nGas  = model.nGas;
            pcoef = model.PoiseuilleCoefficient;
            Tbc   = model.Tbc;
            
            bd = 'Boundary';
            
            qBc = state.(bd).flux;
            pBc = state.(bd).pressure;
            
            p = state.pressure;
            p = mapToBc*p;

            state.(bf).bcFluxDefinition = qBc - pcoef*Tbc*(pBc - p);
            
        end

        function state = updateMassConses(model, state)
            
            nGas = model.nGas;
            div  = model.operators.div;
            
            for igas = 1 : nGas

                j   = state.transferRates{igas}
                src = state.bcSources{igas}
                q   = state.fluxes{igas};

                eqs{igas} = div(q) + j - src;
                
            end

            state.MassConses = eqs;
            
        end

        function state = updateSources(model, state)

            bd = 'Boundary';
            nGas = model.nGas;
            
            for igas = 1 : nGas

                qbc = state.(bd).massFluxes;
                ms{igas} = - mapFromBc*qbc;
                
            end

            state.massSources = ms;

        end


        function state = updateBoundaryMolFractions(model, state)

            bd        = 'Boundary';
            nGas      = model.nGas;
            givenMfBc = model.givenBoundaryVolumeFractions;
            
            pBc = state.(bd).pressure;
            p = state.pressure;            
            p = mapToBc*p;
            
            for igas = 1 : nGas

                mf = state.molFractions{igas};
                mf = mapToBc*mf;
                
                mfBc = state.(bd).molFractions{igas};

                ind = pBc - p < 0;

                bcdefs{igas}      = mfBc - mf;
                bcdefs{igas}(ind) = mfBc(ind) - givenMfBc(ind);
 
            end
            
                
        end
        function state = updateFluxes(model, state)
            
            nGas  = model.nGas;
            pcoef = model.PoiseuilleCoefficient;

            p = state.pressure;
            
            for igas = 1 : nGas

                pigas = state.pressures{igas} 
                xigas = pigas./p;

                q = assembleHomogeneousFlux(model, pcoef, p);

                qs{igas} = assembleUpwindFlux(model, q, xigas);
                
            end
            
            state.fluxes = qs;
        end


        function state = updateControl(model, state)

            bd = 'Boundary';

            bc = model.boundarySetup;
            
            pBc = state.(bd).pressure;
            qBc = state.(bd).flux;

            eqs{1} = bc.pressureMap*pBc - bc.pressureValues;
            eqs{2} = bc.fluxMap*qBc - bc.fluxValues;

            state.(bd).control = vertcat(eqs{:});
            
        end
        
    end

    methods(Static)
        
        function state = updatePressure(model, state)

            nGas = model.nGas;
            
            p = 0*state.pressures{1};

            for igas = 1 : nGas
                p = p + state.pressures{igas};
            end

            state.pressure = p;
            
        end

        function state = updateMolFractions(model, state)
            
            nGas = model.nGas;
            
            p = state.pressure;
            
            for igas = 1 : nGas

                pigas = state.pressures{igas} 
                mfs{igas} = pigas./p;

            end

            state.molFractions = mfs;
            
        end

        
    end
end
    
