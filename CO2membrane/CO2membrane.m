classdef CO2membrane < BaseModel
    
    properties

        Boundary
        Feed
        Permeate
        
        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        perms
        thickness

        %% Advanced parameters
        feedPoiseuilleCoefficient
        permeatePoiseuilleCoefficient
        
    end
    
    methods
        
        function model = CO2membrane(inputparams)
            
            model = model@BaseModel();

            % fdnames = {'G'                        , ...
            %            'perms'                    , ...
            %            'thickness'                , ...
            %            'feedPoiseuilleCoefficient', ...
            %            'permeatePoiseuilleCoefficient'};
            
            % model = dispatchParams(model, inputparams, fdnames);

            % model.Boundary = CO2membraneBoundary(inputparams.Boundary);
            model.Boundary = CO2membraneBoundary([]);
            model.Feed     = CO2Layer([]);
            model.Permeate = CO2Layer([]);
            
            model.constants = PhysicalConstants();
            
            model.gasInd.CO2 = 1;
            model.gasInd.O2  = 2;
            model.gasInd.N2  = 3;
            model.gasInd.Ar  = 4;
            model.nGas = 4;

        end

        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            nGas = model.nGas;
            
            varnames = {};
            % feed fluxes [mol/s]
            varnames{end + 1} = VarName({}, 'feedFluxes', nGas);
            % permeate fluxes [mol/s]
            varnames{end + 1} = VarName({}, 'permeateFluxes', nGas);
            % feed pressures [pascal]
            varnames{end + 1} = VarName({}, 'feedPressures', nGas);
            % permeate pressures [pascal]
            varnames{end + 1} = VarName({}, 'permeatePressures', nGas);
            % transfer rates [mol/s/m]
            varnames{end + 1} = VarName({}, 'transferRates', nGas);
            % feed total pressure [pascal]
            varnames{end + 1} = 'feedPressure';
            % permeate total pressure [pascal]
            varnames{end + 1} = 'permeatePressure';
            % feed boundary sources [mol/s]
            varnames{end + 1} = VarName({}, 'feedBcSources', nGas);
            % permeate boundary sources [mol/s]
            varnames{end + 1} = VarName({}, 'permeateBcSources', nGas);

            
            % Feed mass conservation equations
            varnames{end + 1} = VarName({}, 'feedMassConses', nGas);
            % Permeate mass conservation equations
            varnames{end + 1} = VarName({}, 'permeateMassConses', nGas);
            
            % Boundary condition flux definition equation
            varnames{end + 1} = VarName({'Boundary'}, 'feedBcDefinitions', nGas);
            % Boundary condition flux definition equation
            varnames{end + 1} = VarName({'Boundary'}, 'permeateBcDefinitions', nGas);

            % Boundary condition flux definition equation
            varnames{end + 1} = VarName({'Boundary'}, 'feedBcDefinitions', nGas);
            % Boundary condition flux definition equation
            varnames{end + 1} = VarName({'Boundary'}, 'permeateBcDefinitions', nGas);
            
            % Feed Control 
            varnames{end + 1} = VarName({'Boundary'}, 'feedControls', nGas);            
            % Permeate Control 
            varnames{end + 1} = VarName({'Boundary'}, 'permeateControls', nGas);
            
            model = model.registerVarNames(varnames);

            for igas = 1 : nGas
                
                fn = @CO2membrane.updateFeedMassConses;
                inputvarnames = {VarName({}, 'feedFluxes', nGas, igas), VarName({}, 'transferRates', nGas, igas), VarName({}, 'feedBcSources', nGas, igas)};
                outputvarname = VarName({}, 'feedMassConses', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @CO2membrane.updatePermeateMassConses;
                inputvarnames = {VarName({}, 'permeateFluxes', nGas, igas), VarName({}, 'transferRates', nGas, igas), VarName({}, 'permeateBcSources', nGas, igas)};
                outputvarname = VarName({}, 'permeateMassConses', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @CO2membrane.updateTransferRates;
                inputvarnames = {VarName({}, 'permeatePressures', nGas, igas), VarName({}, 'feedPressures', nGas, igas)};
                outputvarname = VarName({}, 'transferRates', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});                

                fn = @CO2membrane.updateFeedFluxes;
                inputvarnames = {'feedPressure', VarName({}, 'feedPressures', nGas, igas)};
                outputvarname = VarName({}, 'feedFluxes', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @CO2membrane.updatePermeateFluxes;
                inputvarnames = {'permeatePressure', VarName({}, 'permeatePressures', nGas, igas)};
                outputvarname = VarName({}, 'permeateFluxes', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @CO2membrane.updateFeedBCdefinitions;
                inputvarnames = {VarName({}, 'feedPressures', nGas, igas), ...
                                 VarName({'Boundary'}, 'feedPressures', nGas, igas), ...
                                 VarName({'Boundary'}, 'feedFluxes', nGas, igas)};
                outputvarname = VarName({'Boundary'}, 'feedBcDefinitions', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @CO2membrane.updatePermeateBCdefinitions;
                inputvarnames = {VarName({}, 'permeatePressures', nGas, igas), ...
                                 VarName({'Boundary'}, 'permeatePressures', nGas, igas), ...
                                 VarName({'Boundary'}, 'permeateFluxes', nGas, igas)};
                outputvarname = VarName({'Boundary'}, 'permeateBcDefinitions', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @CO2membrane.updateFeedSources;
                inputvarnames = {VarName({'Boundary'}, 'feedFluxes', nGas, igas)};
                outputvarname = VarName({}, 'feedBcSources', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @CO2membrane.updatePermeateSources;
                inputvarnames = {VarName({'Boundary'}, 'permeateFluxes', nGas, igas)};
                outputvarname = VarName({}, 'permeateBcSources', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @CO2membrane.updateFeedControls;
                inputvarnames = {VarName({'Boundary'}, 'feedFluxes', nGas, igas), ...
                                 VarName({'Boundary'}, 'feedPressures', nGas, igas)};
                outputvarname = VarName({'Boundary'}, 'feedControls', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @CO2membrane.updatePermeateControls;
                inputvarnames = {VarName({'Boundary'}, 'permeateFluxes', nGas, igas), ...
                                 VarName({'Boundary'}, 'permeatePressures', nGas, igas)};
                outputvarname = VarName({'Boundary'}, 'permeateControls', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            end

            fn = @CO2membrane.updateFeedPressure;
            inputvarnames = {VarName({}, 'feedPressures', nGas)};
            outputvarname = 'feedPressure';
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            fn = @CO2membrane.updatePermeatePressure;
            inputvarnames = {VarName({}, 'permeatePressures', nGas)};
            outputvarname = 'permeatePressure';
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            varnames = {};
            varnames{end + 1} = {'Boundary', 'feedPressure'};
            varnames{end + 1} = {'Boundary', 'permeatePressure'};

            model = model.setAsExtraVarNames(varnames);
            
        end

        function state = updateFeedMassConses(model, state)
            
            nGas = model.nGas;
            div  = model.operators.div;
            
            for igas = 1 : nGas

                j   = state.transferRates{igas}
                src = state.feedBcSources{igas}
                q   = state.feedFluxes{igas};

                eqs{igas} = div(q) + j - src;
                
            end

            state.feedMassConses = eqs;
            
        end
        
        function state = updatePermeateMassConses(model, state)

            nGas = model.nGas;
            div  = model.operators.div;
            
            for igas = 1 : nGas

                j   = state.transferRates{igas}
                src = state.permeateBcSources{igas}
                q   = state.permeateFluxes{igas};

                eqs{igas} = div(q) - j - src;
                
            end

            state.permeateMassConses = eqs;
            
        end
        
        function state = updateTransferRates(model, state)

            nGas      = model.nGas;
            perms     = model.perms;
            thickness = model.thickness;
            
            for igas = 1 : nGas

                fp = state.feedPressures{igas};
                pp = state.permeatePressures{igas};

                perm = perms(igas);

                js{igas} = perm/thickness*(fp - pp);
                
            end
            
            state.transferRates = js;
        end
        
        function state = updateFeedFluxes(model, state)
            
            nGas  = model.nGas;
            pcoef = model.feedPoiseuilleCoefficient;

            p = state.feedPressure;
            
            for igas = 1 : nGas

                pigas = state.feedPressures{igas} 
                xigas = pigas./p;

                q = assembleHomogeneousFlux(model, pcoef, p);

                qs{igas} = assembleUpwindFlux(model, q, xigas);
                
            end
            
            state.feedFluxes = qs;
        end
        
        function state = updatePermeateFluxes(model, state)
            
            nGas  = model.nGas;
            perms = model.perms;
            pcoef = model.permeatePoiseuilleCoefficient;

            p = state.permeatePressure;
            
            for igas = 1 : nGas

                pigas = state.permeatePressures{igas} 
                xigas = pigas./p;

                q = assembleHomogeneousFlux(model, pcoef, p);

                qs{igas} = assembleUpwindFlux(model, q, xigas);
                
            end
            
            state.permeateFluxes = qs;
            
        end
        
        function state = updateFeedPressure(model, state)

            nGas = model.nGas;
            
            p = 0*state.feedPressures{1};

            for igas = 1 : nGas
                p = p + state.feedPressures{igas};
            end

            state.feedPressure = p;
            
        end
        
        function state = updatePermeatePressure(model, state)

            nGas = model.nGas;
            
            p = 0*state.permeatePressures{1};

            for igas = 1 : nGas
                p = p + state.permeatePressures{igas};
            end

            state.permeatePressure = p;
            
        end

        function state = updateFeedBCdefinitions(model, state)

            bd = 'Boundary';
            nGas = model.nGas;

            for igas = 1 : nGas
                p   = state.feedPressures{igas};
                pBc = state.(bd).feedPressures{igas};
                qBc = state.(bd).feedFluxes{igas};

                p = mapToBc*p;
                
                ind = p - pBc > 0;

                
                
                bcdefs{igas} = ;

            end

            state.(bd).feedBcDefinitions = bcdefs;


        end
        
        function state = updatePermeateBCdefinitions(model, state)

            
        end
        
        function state = updateFeedSources(model, state)

        end
        
        function state = updatePermeateSources(model, state)

        end
        
        function state = updateFeedControls(model, state)

        end
        
        function state = updatePermeateControls(model, state)

        end
        
    end
    
end

