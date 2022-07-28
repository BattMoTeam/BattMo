function schedule = battmoControl2schedule(u, schedule, scaling)
    %% Convert control vector u to schedule 
    
    nc = numel(schedule.control);
    [umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
    c = 0;
    for cs = 1 : nc
        c = c + 1;
        cutoffV = 3.0;
        tup = 0.1;
        Imax =  u(c) *(umax-umin)+umin;
        %schedule.control(cs).Imax = Imax;
        schedule.control(cs).src =@(time,I,E) rampupSwitchControl(time,tup,I,E, Imax, cutoffV);        
    end
    
end