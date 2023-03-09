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

        function [funcstr, statestr, nmstr] = setupCallStringElements(propfunction)
            
            fn = propfunction.fn;
            nmstr = propfunction.modelnamespace;
            nmstr = join(nmstr, '.');
            if ~isempty(nmstr)
                nmstr = nmstr{1};
                statestr = sprintf('state.%s', nmstr);
            else
                statestr = 'state';
            end
            funcstr = func2str(fn);
            funcstr = regexp(funcstr, "\.([^.]*)$", 'tokens');
            funcstr = funcstr{1}{1};
            funcstr = horzcat(nmstr, {funcstr});
            funcstr = join(funcstr, '.');
            funcstr = funcstr{1};
            funcstr = sprintf('model.%s', funcstr);
            
        end

    end
    
    methods(Static)

        function callstr = stdFuncCallSetupFn(propfunction)

            [funcstr, statestr, nmstr] = propfunction.setupCallStringElements();
            
            callstr = sprintf('%s = %s(%s);', statestr, funcstr, statestr);

        end

        function callstr = accumFuncCallSetupFn(propfunction)

            [funcstr, statestr, nmstr] = propfunction.setupCallStringElements();

            state0str = sprintf('state0.%s', nmstr);

            callstr = sprintf('%s = %s(%s, %s, dt);', statestr, funcstr, statestr, state0str);
            
        end

        function callstr = drivingForceFuncCallSetupFn(propfunction)

            [funcstr, statestr, nmstr] = propfunction.setupCallStringElements();
            
            callstr = sprintf('%s = %s(%s, drivingForces);', statestr, funcstr, statestr);
            
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
