classdef PropFunction
    
    properties
        
        varname
        
        inputvarnames
        modelnamespace 
        
        fn % function handler

        functionCallSetupFn %  function handler to setup the function call
    end
    
    methods
        
        function propfunction = PropFunction(varname, fn, inputvarnames, modelnamespace, functionCallSetupFn)
            
            propfunction.varname = varname;
            propfunction.fn = fn;
            propfunction.inputvarnames = inputvarnames;
            propfunction.modelnamespace = modelnamespace;

            if nargin > 4 && ~isempty(functionCallSetupFn)
                propfunction.functionCallSetupFn = functionCallSetupFn;
            else
                propfunction.functionCallSetupFn = @(pfunc) PropFunction.stdFuncCallSetupFn(pfunc);
            end
            
        end
        
        
        function propfunctions = resolveIndex(propfunction)
        % resolve the index for varname (when propfunction.varname is a cell, duplicate the propfunction for each cell entry)
            varname = propfunction.varname;
            varnames = varname.resolveIndex();
            propfunctions = {};
            for ind = 1 : numel(varnames)
                propfunctions{ind} = propfunction;
                propfunctions{ind}.varname = varnames{ind};
            end
        end

        function str = getFunctionCallString(propfunction)
        % We store this function in a function handler that can be changed more dynamically compared to setup a child class

            fn = propfunction.functionCallSetupFn;
            str = fn(propfunction);
            
        end

    end
    
    methods(Static)

        function str = stdFuncCallSetupFn(propfunction)

            fn = propfunction.fn;
            mn = propfunction.modelnamespace;
            mn = join(mn, '.');
            if ~isempty(mn)
                mn = mn{1};
                statename = sprintf('state.%s', mn);
            else
                statename = 'state';
            end
            fnname = func2str(fn);
            fnname = regexp(fnname, "\.(.*)", 'tokens');
            fnname = fnname{1}{1};
            fnname = horzcat(mn, {fnname});
            fnname = join(fnname, '.');
            fnname = fnname{1};

            str = sprintf('%s = model.%s(%s);', statename, fnname, statename);

        end
        

        function str = accumFuncCallSetupFn(propfunction)

            fn = propfunction.fn;
            mn = propfunction.modelnamespace;
            mn = join(mn, '.');
            if ~isempty(mn)
                mn = mn{1};
                statename = sprintf('state.%s', mn);
                statename0 = sprintf('state0.%s', mn);
            else
                statename = 'state';
                statename0 = 'state0';
            end
            fnname = func2str(fn);
            fnname = regexp(fnname, "\.(.*)", 'tokens');
            fnname = fnname{1}{1};
            fnname = horzcat(mn, {fnname});
            fnname = join(fnname, '.');
            fnname = fnname{1};

            str = sprintf('%s = model.%s(%s, %s, dt);', statename, fnname, statename, statename0);

        end
        

    end
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
