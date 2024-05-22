classdef CO2membrane < BaseModel
    
    properties

        Feed
        Permeate
        
        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        permeabilities % structure with permability for each gas
        thickness

        permeabilityValues % vector computed from permeabilities

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

            model.permeabilityValues = permValues;
            
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

        function initstate = setupInitialState(model)

            initstate.Feed     = model.Feed.setupInitialState();
            initstate.Permeate = model.Permeate.setupInitialState();
            
        end


        
        function state = updateSideMassConses(model, state)


            function lstate = update(lstate, side)
            % generic function to update mass cons term for a given side
                nGas = model.(side).nGas;
                G    = model.(side).G; 
                
                for igas = 1 : nGas

                    j   = lstate.transferRates{igas};
                    src = lstate.(side).bcSources{igas};
                    q   = lstate.(side).fluxes{igas};

                    eqs{igas} = G.getDiv(q) - src;

                    switch side
                      case 'Feed'
                        eqs{igas} = eqs{igas} + j;
                      case 'Permeate'
                        eqs{igas} = eqs{igas} - j;
                      otherwise
                        error('not recognized');
                    end
                    
                end

                lstate.(side).massConses = eqs;
            end

            state = update(state, 'Feed');
            state = update(state, 'Permeate');

        end
        

        function state = updateTransferRates(model, state)

            nGas      = model.nGas;
            perms     = model.permeabilityValues;
            thickness = model.thickness;
            
            for igas = 1 : nGas

                fp = state.Feed.pressures{igas};
                pp = state.Permeate.pressures{igas};

                perm = perms(igas);

                js{igas} = perm/thickness*(fp - pp);
                
            end
            
            state.transferRates = js;
        end

        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
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

