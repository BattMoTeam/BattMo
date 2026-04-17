classdef EquivalentCircuitModel < BaseModel

    properties

        nominalcellcapacity
        
        OCP
        R0
        R1
        C1
        R2
        C2
        
        initSOC
        initOverpotential2
        lowerVoltageCutoff

        I
        totalTime
        
        %% helpers
        OCPfunc
        Ifunc 
    end

    methods

        function model = EquivalentCircuitModel(inputparams)

            model = model@BaseModel();

            %% Setup the model using the input parameters
            
            fdnames = {'nominalcellcapacity', ...
                       'OCP'                , ...
                       'R0'                 , ...
                       'R1'                 , ...
                       'C1'                 , ...
                       'R2'                 , ...
                       'C2'                 , ...
                       'initSOC'            , ...
                       'initOverpotential2' , ...
                       'lowerVoltageCutoff' , ...
                       'I'                  , ...
                       'totalTime'};

            model = dispatchParams(model, inputparams, fdnames);

            model.OCPfunc = setupFunction(model.OCP);
            model.Ifunc   = setupFunction(model.I);
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % Current
            varnames{end + 1} = 'I';
            %
            varnames{end + 1} = 'OCP';
            %
            varnames{end + 1} = 'U';
            
            varnames{end + 1} = 'UR';            
            varnames{end + 1} = 'U1';
            varnames{end + 1} = 'U2';
            varnames{end + 1} = 'SOC';

            varnames{end + 1} = 'dU1dt';            
            varnames{end + 1} = 'dU2dt';
            varnames{end + 1} = 'dSOCdt';

            model = model.registerVarNames(varnames);

            fn = @EquivalentCircuitModel.updateUR;
            inputnames = {'I'};
            model = model.registerPropFunction({'UR', fn, inputnames});

            fn = @EquivalentCircuitModel.updatedU1dt;
            inputnames = {'U1', 'I'};
            model = model.registerPropFunction({'dU1dt', fn, inputnames});
                        
            fn = @EquivalentCircuitModel.updatedU2dt;
            inputnames = {'U2', 'I'};
            model = model.registerPropFunction({'dU2dt', fn, inputnames});

            fn = @EquivalentCircuitModel.updatedSOCdt;
            inputnames = {'I'};
            model = model.registerPropFunction({'dSOCdt', fn, inputnames});

            fn = @EquivalentCircuitModel.updateU;            
            inputnames = {'SOC', 'U1', 'U2', 'UR'};
            model = model.registerPropFunction({'U', fn, inputnames});            

        end

        function state = updateUR(model, state)

            state.UR = model.R0*state.I;
            
        end
        
        function state = updatedU1dt(model, state)
            
            state.dU2dt = model.updatedUdt(state.U2, model.R2, model.C2, state.I);
            
        end
        
        function state = updatedU2dt(model, state)

            state.dU1dt = model.updatedUdt(state.U1, model.R1, model.C1, state.I);
            
        end
        
        function state = updatedSOCdt(model, state)

            Qmax = model.nominalcellcapacity;
            
            state.dSOCdt = -(1 / Qmax) * state.I;
            
        end
        
        function state = updateU(model, state)

            OCP = model.OCPfunc(state.SOC);
            state.U = OCP - (state.U1 + state.U2 + state.UR);
            
        end

        function state0 = setupInitialCondition(model, SOC0)

            state0.U1  = 0;
            state0.U2  = 0;
            state0.SOC = SOC0;
            
        end

        function state = setupStateFromY(model, y)

            state.U1  = y(1);
            state.U2  = y(2);
            state.SOC = y(3);
            
        end

        function y = setupYfromState(model, state)

            y(1) = state.U1;
            y(2) = state.U2;
            y(3) = state.SOC;
            
        end
        
        function y = setupFfromState(model, state)

            y(1) = state.dU1dt;
            y(2) = state.dU2dt;
            y(3) = state.dSOCdt;
            y = y';
            
        end
        
        function f = ode(model, t, y)

            state = model.setupStateFromY(y);

            state.I = model.Ifunc(t);
            
            state = model.evalVarName(state, 'dU1dt');
            state = model.evalVarName(state, 'dU2dt');
            state = model.evalVarName(state, 'dSOCdt');
            
            f = model.setupFfromState(state);
            
        end
        
        function [t, U, I] = solve(model)

            if isempty(model.computationalGraph)
                model = setupComputationalGraph(model);
            end

            % setup initial condition
            state0 = model.setupInitialCondition(model.initSOC);
            y0 = model.setupYfromState(state0);

            % setup time span
            tspan = [0, model.totalTime];

            % solve ode
            [t_out, y_out] = ode45(@(t, y) model.ode(t, y), tspan, y0);

            for i = 1:length(t_out)
                
                state = model.setupStateFromY(y_out(i, :));
                I(i) = model.Ifunc(t_out(i));
                state.I = I(i);
                state = model.evalVarName(state, 'U');
                U(i) = state.U;
            end

            t = t_out;
        end
        
    end

    methods (Static)

        function dUdt = updatedUdt(U, R, C, I)
            dUdt = -(1 / (R * C)) * U + (1 / C) * I;
        end
        
    end
    
end
