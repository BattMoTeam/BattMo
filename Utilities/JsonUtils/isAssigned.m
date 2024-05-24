function status = isAssigned(value)

    if isa(value, 'UnAssigned')
        status = false;
    else
        status = true;
    end
    
end
