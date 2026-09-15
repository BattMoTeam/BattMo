function [fn_handler, fn] = setupFunction(jsonstruct)

    switch jsonstruct.functionFormat
        
      case 'tabulated'

        numberOfArguments = numel(jsonstruct.argumentList);
        
        switch numberOfArguments

          case 1
            
            fn = TabulatedFunction1D(jsonstruct);
            
          case 2
            
            error('bilinear interpolate function not supported yet');
            
          otherwise

            error('the given number of arguments %d is not supported', numberOfArguments);
        end

        fn_handler = @(x) fn.eval(x);

        return
        
      case 'string expression'

        fn = FormulaFunction(jsonstruct);
        fn_handler = @(varargin) fn.eval(varargin{:});

        return
        
      case 'named function'
        
        fn = NamedFunction(jsonstruct);
        fn_handler = @(varargin) fn.eval(varargin{:});

        return
        
      case 'constant'

        fn = ConstantFunction(jsonstruct);
        fn_handler = @(varargin) fn.eval(varargin{:});
        
      case 'csv'

        data = readmatrix(jsonstruct.filename);
        
        argumentList = jsonstruct.argumentList;

        jsonstruct = struct('functionFormat', 'tabulated' , ...
                            'dataX'         , data(: , 1) , ...
                            'dataY'         , data(: , 2));

        jsonstruct.argumentList = argumentList; 
        
        [fn_handler, fn] = setupFunction(jsonstruct);

        return
        
      otherwise

        error('function format %s not recognized', jsonstruct.functionFormat);

    end
        
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
