classdef CO2membrane < BaseModel
    
    properties

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
            model.Feed     = CO2membraneSide([]);
            model.Permeate = CO2membraneSide([]);
            
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
        

        

    end
    
end

