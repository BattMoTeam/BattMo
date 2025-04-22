function [fn, fn_handler] = setupFunction(jsonstruct)

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

      otherwise

        error('function format not recognized');
        
    end
        
end
