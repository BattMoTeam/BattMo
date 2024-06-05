classdef StringExpFunction
% The function is given by a string
    
    properties
        
        stringExpression % The name of the variable should match those given in the argument list 

        %% helper

        formattedStringExpression 
        
    end

    methods

        function fn = Function(jsonstruct)

            fdnames = {'stringExpression'};

            fn = dispatchParams(fn, jsonstruct, fdnames);

            nArgs = fn.numberOfArguments;
            args  = fn.argumentList;

            fmtstr = '';

            for iarg = 1 : nArgs

                fmtstr = sprintf('%s; %s = varargin{%d};', fmtstr, args{iarg}, iarg);
                
            end

            fmtstr = sprintf('%s; y = %s;', fmtstr, fn.stringExpression);

            fn.formattedStringExpression = fmtstr;
            
        end
        
        function y = eval(fn, varargin)

            eval(fn.formattedStringExpression);
            
        end

        
    end
    
    
end
