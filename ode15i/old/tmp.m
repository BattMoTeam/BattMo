   y = getOdeVariables(state);
   state = stateFromOdeVarables(y);


function odeFun(t, y, dy)
    %if(derivget)
        y = initVariablesAD(y);
    %end
    state = getOdeVariableToState(y);
    state0 = state
    dt = 1e100 % dummy
    %% all transport terms
    eq1 = model.getEquations(state,state,dt,force, []);
    %if(derivget)
    dy = initVariablesAD(dy);
    %end
    eq2_tmp = model.getAccumulationterm(state);%% implemet
    eq2 = eq2*value(dy);
    eq = value(eq1) + eq2;
    %if(derivget)
        f = eq
        dfdy = eq1.jac
        dfdydot = eq2_tmp.jac;
    %else
    %    f = eq;
    %end
end
    
