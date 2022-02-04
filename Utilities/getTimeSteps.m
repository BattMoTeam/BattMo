function times = getTimeSteps(dt0, n, total, fac)
    dt=[];
    dt1 = total/n;
    dt0_org = dt0;
    nstart= ceil(-log(dt0_org)/log(fac));
    dt = [dt; repmat(dt1, nstart, 1).*fac.^[1:nstart]'/fac^nstart];
    % Time scaling can be adding using variable tfac
    dt = [dt; repmat(dt1, n, 1)]; 
    times = [0; cumsum(dt)];
end