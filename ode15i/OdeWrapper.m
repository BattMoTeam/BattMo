classdef OdeWrapper
    properties
        model
        schedule
        initstate
        forces
        VariableNames
        VariableSizes
        initVariables
        currentState
    end
    methods
        function wrap = OdeWrapper(model,initstate,forces,schedule)
            %Constructor. We initialize the ode variables needed by ode15i
            %using a Battery model, an initial state and precalculated
            %forces
            wrap.model=model;
            wrap.initstate=initstate;
            wrap.schedule=schedule;
            %Assumes forces to be a cell containing src and IEswitch names
            %+ values
            wrap.forces=cell2struct({forces{2},forces{4}}, {forces{1},forces{3}},2);
            wrap.VariableNames=model.primaryVariableNames;

            %% Setup initial state variables
            %NB! We assume dy/dt at t=0 to be 0 !!

            function [val, valSize] = dummy(x)
                val=model.getProp(initstate,x);
                valSize=length(val);
            end
            [inity,dims] = cellfun(@dummy , wrap.VariableNames , 'UniformOutput' , false);
            wrap.VariableSizes=cell2mat(dims);

            inity = cell2mat(inity); % ode15i requires vector format
            initdy = inity.* 0; %Should be a better way to do this
            wrap.initVariables = struct('y',inity,'dy',initdy);
            wrap.currentState=initstate;
        end

        function [eq,other,accum]=odeEqs(wrap,t,y,dy)
            %Generate the governing equations of the problem and set them
            %up to the compatible with the inputs of ode15i

            %% Initialize AD variables
            y = initVariablesADI(y);
            dy = initVariablesADI(dy);

            state=wrap.getStateFromVariables(y); 

            dt = 1; % dummy
            
            %% all transport terms 
            state.time = t;
            [other_tmp,state]= wrap.model.getEquations(state,state,dt,wrap.forces, 'ResOnly',true); %X
            other=EquationVector(other_tmp.equations);
          
            accum_tmp = wrap.model.getAccumTerms(state);
            accum_tmp = EquationVector(accum_tmp);
            accum=accum_tmp.jac{1}*dy;
            eq = value(accum) + value(other);

            wrap.currentState=state;
        end

        function state=getStateFromVariables(wrap,y)
            %Generate a state from a vector of ode variables, these should
            %correspond to the primary variables of the model.

            ind=1;
            %Use current state to copy all variables that are not among the
            %primary variables
            state=wrap.currentState;
            for i=1:length(wrap.VariableSizes)
                state=wrap.model.setNewProp(state, ...
                    wrap.VariableNames{i}, ...
                    y(ind:ind + wrap.VariableSizes(i)-1));
                ind = ind + wrap.VariableSizes(i);
            end
        end

        function [res, ret_state] = solve(wrap, varargin)
            %% Reset
            wrap.currentState=wrap.initstate;
            %Index of voltage
            E_ind=sum(wrap.VariableSizes(1:5));

            %% Functions used in ode15i
            function f = odeFun(t,y,dy)
                [f,~,~] = wrap.odeEqs(t,y,dy);
            end

            function [dfdy , dfdyp] = jacFun(t,y,dy)
                [~,tmp_dfdy,tmp_dfdyp] = wrap.odeEqs(t,y,dy);
                dfdy=tmp_dfdy.jac{1};
                dfdyp=tmp_dfdyp.jac{1};
            end
            
            %% Stopping condition

            function [val,isterminal,direction]=cutoffEvent(t,y,varargin)
                val = (y(E_ind) >= wrap.model.Control.lowerCutoffVoltage + 1e-4) & (y(E_ind) <= wrap.model.Control.upperCutoffVoltage + 1);
                isterminal = 1;
                direction = 0;
            end

            %% Build odeset input options
            options = odeset('Jacobian',@jacFun,'Events',@cutoffEvent, "Stats", "on");
            %options=odeset();
            %Solve problem
            res = ode15i(@odeFun, cumsum(wrap.schedule.step.val),wrap.initVariables.y,wrap.initVariables.dy,options);
            ret_state= wrap.getStateFromVariables(res.y);
            ret_state.E=res.y(E_ind,:);
        end
    end
end