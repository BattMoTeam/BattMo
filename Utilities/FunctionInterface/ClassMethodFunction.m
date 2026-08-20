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
            
            if ~isempty(fn.parameters)
                
                fn.obj = feval(fn.className, fn.parameters);
                
            else
                
                mc = meta.class.fromName(fn.className);
                idx = find(strcmp({mc.MethodList.Name}, fn.methodName), 1);
                isStatic = ~isempty(idx) && mc.MethodList(idx).Static;

                if ~isStatic
                    fn.obj = feval(fn.className);
                else
                    fn.obj = [];
                end
                
            end

            if ~isempty(fn.obj)
                fn.functionHandler = @(varargin) fn.obj.(fn.methodName)(varargin{:});
            else
                % static/class method handle like @() ClassName.methodName(...)
                fn.functionHandler = str2func([fn.className '.' fn.methodName]);
            end
            
        end
        
        function y = eval(fn, varargin)

            y = fn.functionHandler(varargin{:});

        end
        
    end
    
    
end
