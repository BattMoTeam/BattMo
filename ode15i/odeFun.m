function f_ret=odeFun(t, y, dy,state0,dims,model,forces)
    [f,~,~]=odeEqs(t,y,dy,state0,dims,model,forces);
    f_ret=f.val;
end