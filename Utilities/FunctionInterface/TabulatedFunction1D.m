classdef TabulatedFunction1D < TabulatedFunction
% One-argument tabulated function. The value is obtained using linear interpolation.

    properties

        dataX % Array for input values
        
        dataY % Array for output values (same size as dataX)
        
    end

    methods

        function fn = TabulatedFunction1D(jsonstruct)

            fn = fn@TabulatedFunction(jsonstruct);

            fdnames = {'dataX', ...
                       'dataY'};

            fn = dispatchParams(fn, jsonstruct, fdnames);

        end

        function y = eval(fn, x)

            y = interpTable(fn.dataX, fn.dataY, x);
            
        end

    end
    

end
