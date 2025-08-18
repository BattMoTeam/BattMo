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
            
            error('the given number of arguments is not supported');
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

      case 'csv'

        data = readmatrix(jsonstruct.filename);

        jsonstruct = struct('functionFormat', 'tabulated'            , ...
                            'argumentList'  , jsonstruct.argumentList, ...
                            'dataX'         , data(:                 , 1), ...
                            'dataY'         , data(:                 , 2));
        
        [fn_handler, fn] = setupFunction(jsonstruct);

        return
        
      otherwise

        error('function format not recognized');
        
    end
        
end
