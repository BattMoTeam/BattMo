function [eq, eq1,eq2_tmp]=odeEqs(t,y,dy,state0,dims,model,forces)
    y = initVariablesADI(y);
    %Look at problem output from getEquations. Identify members of list
    %in state and associate each one with an element in the array y.
    %The y -> state function is then simply the inverse mapping
    %Question: What to do with remaining variables?
    %Solved X
    state = getStateFromVariables(y,model,dims,state0); 
    %For the moment
    state0 = state;
    dt = 1; % dummy
    
    %% all transport terms

    %Use default forces
    %forces=model.getValidDrivingForces();

    %Reshape

    %Temp, I don't know what to do with this. IESwitch gives error
    %model.Control.controlPolicy='None';
    tmp = model.getEquations(state,state,dt,forces, 'ResOnly',true); %X
    eq1=EquationVector(tmp.equations);

    dy = initVariablesADI(dy); %X

    %Check problem object returned by getEquations
    %Charge, EI, Control do not accumulate!
    %Look at AccumTerm for mass + solid diffusion
    %These can become differentiable through AD. This is what we evaluate
    %Solved X
    tmp = model.getAccumTerms(state);
    eq2_tmp=EquationVector(tmp);
    
    eq2 = eq2_tmp.*value(dy);
    eq = value(eq1) + eq2; %Value: See AD object
end