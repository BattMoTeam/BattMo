
function [dfdy,dfdyp]=JACodeFun(t,y,dy,state0,dims,model,forces)
    [~,eq1,eq2]=odeEqs(t,y,dy,state0,dims,model,forces);
    dfdy=eq1.jac{1};
    dfdyp=eq2.jac{1};
end   