function jsonstruct = battMojsondecode(jsonstring)
%% wrapper over jsondecode to avoid compatibility issues between Matlab and Octave

    if mrstPlatform('octave')

        jsonstruct = jsondecode(jsonstring, 'makeValidName', false);
        
    else

        jsonstruct = jsondecode(jsonstring);
        
    end

end
