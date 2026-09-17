classdef FormulaFunction < Function
% The function is given by a string
    
    properties
        
        formula       % Expression that should be evaluated
        variableNames % Names of the variables in the expression formula, which corresponds to the order given in the argument list

        %%
        
        outputVarName = 'y' % This is the string that will be used to store the value. By default it is 'y' but we have the possibility to change it. NOTE : As the implementation is done now, it should not conflict with other variables.
        
        %% helper

        formattedStringExpression
        
    end

    methods

        function fn = FormulaFunction(jsonstruct)

            fn = fn@Function(jsonstruct);

            if isfield(jsonstruct, 'expression')

                jsonstruct = setDefaultStructField(jsonstruct, {'expression', 'language'}, 'matlab');
                jsonstruct = setDefaultStructField(jsonstruct, {'expression', 'variableNames'}, jsonstruct.argumentList);
                
                assert(isEqualStructField(jsonstruct, {'expression', 'language'}, 'matlab'), 'The language should be matlab');

                fn.formula       = jsonstruct.expression.formula;
                fn.variableNames = jsonstruct.expression.variableNames;
                
            elseif isfield(jsonstruct, 'expressions')

                nexpr = numel(jsonstruct.expressions);
                iexpr = 1;

                found = false;
                
                while ~found && iexpr <= nexpr
                
                    % depending on json parsing, we get an struct-array or a cell-array 
                    if iscell(jsonstruct.expressions)
                        expr = jsonstruct.expressions{iexpr};
                    else
                        expr = jsonstruct.expressions(iexpr);
                    end

                    if strcmp(expr.language, 'matlab')

                        found = true;
                        expr = setDefaultStructField(expr, {'variableNames'}, jsonstruct.argumentList);

                        fn.formula       = expr.formula;
                        fn.variableNames = expr.variableNames;
                        
                    end

                    iexpr = iexpr + 1;
                    
                end

                assert(found, 'could not find an expression for matlab');
                
            else
                
                error('missing expressions')
                
            end
            
            varnames = fn.variableNames;
            nvars    = numel(varnames);

            assert(nvars == numel(fn.argumentList), 'The number of variables should be the same as the number of argument in the argument list')
            
            fmtstrs  = {};

            for ivar = 1 : nvars

                fmtstrs{ivar} = sprintf('%s = varargin{%d}', varnames{ivar}, ivar);
                
            end

            fmtstr = join(fmtstrs);

            fmtstr = sprintf('%s; %s = %s;', fmtstr{1}, fn.outputVarName, fn.formula);

            fn.formattedStringExpression = fmtstr;
            
        end
        
        function y = eval(fn, varargin)

            eval(fn.formattedStringExpression);
            
        end

        
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
