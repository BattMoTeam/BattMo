classdef StringExpFunction < Function
% The function is given by a string
    
    properties
        
        stringExpression % Expression that should be evaluated
        variableNames % Names of the variables in the expression stringExpression, which corresponds to the order given in the argument list

        %%
        
        outputVarName = 'y' % This is the string that will be used to store the value. By default it is 'y' but we have the possibility to change it. NOTE : As the implementation is done now, it should not conflict with other variables.
        
        %% helper

        formattedStringExpression 
    end

    methods

        function fn = StringExpFunction(jsonstruct)

            fn = fn@Function(jsonstruct);
            
            fdnames = {'stringExpression', ...
                       'variableNames'};

            fn = dispatchParams(fn, jsonstruct, fdnames);

            varnames = fn.variableNames;
            nvars    = numel(varnames);

            assert(nvars == numel(fn.argumentList), 'The number of variables should be the same as the number of argument in the argument list')
            
            fmtstrs  = {};

            for ivar = 1 : nvars

                fmtstrs{ivar} = sprintf('%s = varargin{%d}', varnames{ivar}, ivar);
                
            end

            fmtstr = join(fmtstrs);

            fmtstr = sprintf('%s; %s = %s;', fmtstr{1}, fn.outputVarName, fn.stringExpression);

            fn.formattedStringExpression = fmtstr;
            
        end
        
        function y = eval(fn, varargin)

            eval(fn.formattedStringExpression);
            
        end

        
    end
    
    
end
