classdef CO2membrane < BaseModel
    
    properties

        Feed
        Permeate
        
        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        permeabilities % structure with permability
        thickness

        permValues % vector computed from permeabilities

    end
    
    methods
        
        function model = CO2membrane(inputparams)
            
            model = model@BaseModel();
            
            fdnames = {'permeabilities', ...
                       'thickness'};

            model = dispatchParams(model, inputparams, fdnames);
            
            model.Feed     = CO2membraneSide(inputparams.Feed);
            model.Permeate = CO2membraneSide(inputparams.Permeate);
            
            model.constants = PhysicalConstants();

            model = CO2membrane.setupGasStructures(model);

            nGas   = model.nGas;
            gasInd = model.gasInd;
            
            fdnames = fieldnames(model.permeabilities);
            assert(numel(fdnames) == nGas, 'problem in setup');
            
            for igas = 1 : nGas
                fdname = fdnames{igas};
                permValues(gasInd.(fdname)) = model.permeabilities.(fdname);
            end

            model.permValues = permValues;
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            nGas = model.nGas;
            
            varnames = {};
            % transfer rates [mol/s/m]
            varnames{end + 1} = VarName({}, 'transferRates', nGas);
            
            model = model.registerVarNames(varnames);

            sides = {'Feed', 'Permeate'};

            for igas = 1 : nGas
                
                for  iside = 1 : numel(sides)

                    side = sides{iside};
                    
                    fn = @CO2membrane.updateSideMassConses;
                    inputvarnames = {VarName({side}, 'fluxes', nGas, igas)   , ...
                                     VarName({}, 'transferRates', nGas, igas), ...
                                     VarName({side}, 'bcSources', nGas, igas)};
                    outputvarname = VarName({side}, 'massConses', nGas, igas);
                    model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                end
                
                fn = @CO2membrane.updateTransferRates;
                inputvarnames = {VarName({'Feed'}, 'pressures', nGas, igas), ...
                                 VarName({'Permeate'}, 'pressures', nGas, igas)};
                outputvarname = VarName({}, 'transferRates', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});                

                
            end

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
        

        function state = updateTransferRates(model, state)

            nGas      = model.nGas;
            perms     = model.permeabilities;
            thickness = model.thickness;
            
            for igas = 1 : nGas

                fp = state.feedPressures{igas};
                pp = state.permeatePressures{igas};

                perm = perms(igas);

                js{igas} = perm/thickness*(fp - pp);
                
            end
            
            state.transferRates = js;
        end
        
    end

    methods (Static)

        function model = setupGasStructures(model)

            gasInd.CO2 = 1;
            gasInd.O2  = 2;
            gasInd.N2  = 3;
            gasInd.H2O = 4;
            gasInd.Ar  = 5;

            nGas = numel(fieldnames(gasInd));

            model.gasInd = gasInd;
            model.nGas   = nGas;
            
        end

    end
    
end

