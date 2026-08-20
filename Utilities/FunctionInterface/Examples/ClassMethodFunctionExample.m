classdef ClassMethodFunctionExample

    properties

        param1
        param2
        
    end
    
    methods

        function fn = ClassMethodFunctionExample(parameters)

            if nargin < 1
                parameters = [];
            end
            
            fn = setupParameters(fn, parameters);
            
        end

        function y = method1(fn, x)

            y = fn.param1 + x;
            
        end
        
        function z = method2(fn, x, y)

            z = fn.param2*x + y;
            
        end

        function y = method3(fn, x)

            y = 5*x;
            
        end
        
    end

    methods(Static)

        function y = method4(x)

            y = 2*x;

        end
            

    end

end
