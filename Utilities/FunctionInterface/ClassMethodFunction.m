classdef ClassMethodFunction < Function
% The function is implemented and available in the path. It then called using its name.
    properties

        %%
        %
        className
        methodName
        parameters % list of parameters

        obj
        
        %% helper
        %
        functionHandler
        
    end

    methods


        function fn = ClassMethodFunction(jsonstruct)

            fn = fn@Function(jsonstruct);
            
            fdnames = {'className', ...
                       'methodName', ...
                       'parameters'};

            fn = dispatchParams(fn, jsonstruct, fdnames);
            fn = fn.setupFunctionHandler();
            
        end

        function fn = setupFunctionHandler(fn)
            
            % Create an instance or function handle for the target class/method.
            if isempty(fn.className)
                error('ClassMethodFunction:MissingClassName','No className specified');
            end
            if isempty(fn.methodName)
                error('ClassMethodFunction:MissingMethodName','No methodName specified');
            end
            % Try to instantiate object using supplied parameters (if any).
            % try
                if ~isempty(fn.parameters)
                    fn.obj = feval(fn.className, fn.parameters);
                else
                    fn.obj = feval(fn.className);
                end
            % catch
            %     % If construction fails, assume the target is a static/class method.
            %     fn.obj = [];
            % end
            % Set the function handler to call the method on the object (instance method)
            % or a class/static method.
            if ~isempty(fn.obj)
                fn.functionHandler = @(varargin) fn.obj.(fn.methodName)(varargin{:});
            else
                % static/class method handle like @() ClassName.methodName(...)
                fn.functionHandler = str2func([fn.className '.' fn.methodName]);
            end
        end
        
        function y = eval(fn, varargin)
            % Use the prepared function handler if available.
            if isempty(fn.functionHandler)
                fn = fn.setupFunctionHandler();
            end
            if isa(fn.functionHandler, 'function_handle')
                y = fn.functionHandler(varargin{:});
            else
                % Fallback to calling a method on this object (legacy behavior).
                y = fn.(fn.methodName)(varargin{:});
            end
        end
        
    end
    
    
end
