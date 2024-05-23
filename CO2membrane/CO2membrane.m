classdef CO2membrane < BaseModel
    
    properties

        Feed
        Permeate
        
        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        permeabilities % structure with permability for each gas
        thickness

    end
    
    properties (SetAccess = immutable)

        %% helpers
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
                    
                    fn = @CO2membrane.updateMassSources;
                    inputvarnames = {VarName({}, 'transferRates', nGas, igas), ...
                                     VarName({side}, 'bcSources', nGas, igas)};
                    outputvarname = VarName({side}, 'massSources', nGas, igas);
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

            function lstate = update(lstate, side)
            % generic function to update mass source term for a given side
                
                nGas = model.(side).nGas;
                
                for igas = 1 : nGas
                    
                    % transfer rate is positive from feed to permeate
                    jrate = lstate.transferRates{igas};
                    bcsrc = lstate.(side).bcSources{igas};

                    msrcs{igas} = bcsrc;

                    switch side
                      case 'Feed'
                        msrcs{igas} = msrcs{igas} - jrate;
                      case 'Permeate'
                        msrcs{igas} = msrcs{igas} + jrate;
                      otherwise
                        error('not recognized');
                    end
                    
                end

                lstate.(side).massSources = msrcs;
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

        function [state, report] = updateState(model, state, problem, dx, drivingForces)       

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);

            % cap pressures

            bd = 'Boundary';
            sides = {'Feed', 'Permeate'};

            for iside = 1 : numel(sides)

                side = sides{iside};
                
                nGas = model.(side).nGas;
                
                sum_side    = 0;
                sum_bd_side = 0;
                
                for igas = 1 : nGas

                    state.(side).molFractions{igas} = max(0, state.(side).molFractions{igas});
                    state.(side).molFractions{igas} = min(1, state.(side).molFractions{igas});

                    sum_side = sum_side + state.(side).molFractions{igas};
                    
                    state.(side).(bd).molFractions{igas} = max(0, state.(side).(bd).molFractions{igas});
                    state.(side).(bd).molFractions{igas} = min(1, state.(side).(bd).molFractions{igas});
                    
                    sum_bd_side = sum_bd_side + state.(side).(bd).molFractions{igas};

                end
                
                for igas = 1 : nGas
                    state.(side).molFractions{igas}      = state.(side).molFractions{igas}./sum_side;
                    state.(side).(bd).molFractions{igas} = state.(side).(bd).molFractions{igas}./sum_bd_side;
                end                
            end
            
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

