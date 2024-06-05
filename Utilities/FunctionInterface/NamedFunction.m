classdef NamedFunction
% The function is implemented and available in the path. It then called using its name.
    properties
        
        functionName

        %% helper

        functionHandler
        
    end

    methods


        function fn = NamedFunction(jsonstruct)

            fdnames = {'functionName'};

            fn = dispatchParams(fn, jsonstruct, fdnames);

            fn.functionHandler = str2func(fn.functionName);
            
        end
        
        function [varargout] = eval(fn, varargin)

            varargout = cell(1, nargout);
            [varargout{:}] = fn.functionHandler(varargin{:});
            
        end
        
    end
    
    
end
