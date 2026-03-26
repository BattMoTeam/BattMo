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
