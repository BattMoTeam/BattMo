function ret=EquationVector(eq)
    ret=[];
    for i=1:length(eq)
        ret=[ret;eq{i}];
    end
end