classdef CO2capture < BaseModel
    
    properties

        Feed
        Permeate
        
        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        permeances % structure with permability for each gas
        thickness

    end
    
    properties (SetAccess = immutable)

        %% helpers
        permeabilityValues % vector computed from permeabilities

    end
    
    methods
        
        function model = CO2capture(inputparams)
            
            model = model@BaseModel();
            
            fdnames = {'permeabilities', ...
                       'thickness'};

            model = dispatchParams(model, inputparams, fdnames);
            
            model.Feed     = CO2captureChannel(inputparams.Feed);
            model.Permeate = CO2captureChannel(inputparams.Permeate);
            
            model.constants = PhysicalConstants();

            model = CO2capture.setupGasStructures(model);

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
            % transfer fluxes [mol/s/m]
            varnames{end + 1} = VarName({}, 'crossFluxes', nGas);
            
            model = model.registerVarNames(varnames);

            channels = {'Feed', 'Permeate'};

            for igas = 1 : nGas
                
                for  ichannel = 1 : numel(channels)

                    channel = channels{ichannel};
                    
                    fn = @CO2capture.updateMassSources;
                    inputvarnames = {VarName({}, 'crossFluxes', nGas, igas), ...
                                     VarName({channel}, 'bcSources', nGas, igas)};
                    outputvarname = VarName({channel}, 'massSources', nGas, igas);
                    model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                end
                
                fn = @CO2capture.updateTransferRates;
                inputvarnames = {VarName({'Feed'}, 'pressures', nGas, igas), ...
                                 VarName({'Permeate'}, 'pressures', nGas, igas)};
                outputvarname = VarName({}, 'crossFluxes', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});                

                
            end

        end

        function initstate = setupInitialState(model)

            initstate.Feed     = model.Feed.setupInitialState();
            initstate.Permeate = model.Permeate.setupInitialState();
            
        end

        function state = addVariables(model, state)

        % Given a state where only the primary variables are defined, this
        % functions add all the additional variables that are computed in the assembly process and have some physical
        % interpretation.
        %
        % To do so, we use getEquations function and sends dummy variable for state0, dt and drivingForces

            dt            = 1;
            state0        = state;
            drivingForces = model.getValidDrivingForces();

            [~, state] = getEquations(model, state0, state, dt, drivingForces, 'ResOnly', true);

        end

        
        
        function state = updateMassSources(model, state)

            function lstate = update(lstate, channel)
            % generic function to update mass source term for a given channel
                
                nGas = model.(channel).nGas;
                
                for igas = 1 : nGas
                    
                    % transfer rate is positive from feed to permeate
                    jrate = lstate.crossFluxes{igas};
                    bcsrc = lstate.(channel).bcSources{igas};

                    msrcs{igas} = bcsrc;

                    switch channel
                      case 'Feed'
                        msrcs{igas} = msrcs{igas} - jrate;
                      case 'Permeate'
                        msrcs{igas} = msrcs{igas} + jrate;
                      otherwise
                        error('not recognized');
                    end
                    
                end

                lstate.(channel).massSources = msrcs;
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
            
            state.crossFluxes = js;
        end

        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)       

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);

            % cap pressures

            bd = 'Boundary';
            channels = {'Feed', 'Permeate'};

            for ichannel = 1 : numel(channels)

                channel = channels{ichannel};
                
                nGas = model.(channel).nGas;
                
                sum_channel    = 0;
                sum_bd_channel = 0;
                
                for igas = 1 : nGas

                    state.(channel).molFractions{igas} = max(0, state.(channel).molFractions{igas});
                    state.(channel).molFractions{igas} = min(1, state.(channel).molFractions{igas});

                    sum_channel = sum_channel + state.(channel).molFractions{igas};
                    
                    state.(channel).(bd).molFractions{igas} = max(0, state.(channel).(bd).molFractions{igas});
                    state.(channel).(bd).molFractions{igas} = min(1, state.(channel).(bd).molFractions{igas});
                    
                    sum_bd_channel = sum_bd_channel + state.(channel).(bd).molFractions{igas};

                end
                
                for igas = 1 : nGas
                    state.(channel).molFractions{igas}      = state.(channel).molFractions{igas}./sum_channel;
                    state.(channel).(bd).molFractions{igas} = state.(channel).(bd).molFractions{igas}./sum_bd_channel;
                end                
            end
            
        end
        
    end

    methods (Static)

        function model = setupGasStructures(model, gasSpecies)
            
            for igas = 1 : numel(gasSpecies)
                gasInd.(gasSpecies{igas}) = igas;
            end

            nGas = numel(fieldnames(gasInd));

            model.gasInd = gasInd;
            model.nGas   = nGas;
            
        end

    end
    
end

