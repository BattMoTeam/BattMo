classdef Function

% Parent class for an input function. The function can be given using different format. Each format is implemented with
% its own class. Use the separates function setupFunction to initiate the function

    properties
        
        functionFormat % String that can take values
                       % - 'tabulated',
                       % - 'string expression',
                       % - 'named function',
                       % - 'constant'

        argumentList % cell array of string which describes the expected argument of the function. Mainly for documentation at the moment or to infer number of argument

        %% helper

        numberOfArguments
        
    end

    methods


        function fn = Function(jsonstruct)

            fdnames = {'functionFormat', ...
                       'argumentList'};

            fn = dispatchParams(fn, jsonstruct, fdnames);

            fn.numberOfArguments = numel(fn.argumentList);

        end

        function y = eval(fn, varargin)
            
            error('virtual function');
            
        end

    end
    
    
end
