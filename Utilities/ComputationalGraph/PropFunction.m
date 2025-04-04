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
            elseif isempty(fn)
                propfunction.functionCallSetupFn = @(pfunc) [];
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

            fnsetup = propfunction.functionCallSetupFn;
            str = fnsetup(propfunction);
            
        end

        function [funcstr, statestr, nmstr] = setupCallStringElements(propfunction)
            
            fn = propfunction.fn;
            nmstr = propfunction.modelnamespace;
            nmstr = strjoin(nmstr, '.');
            if ~isempty(nmstr)
                statestr = sprintf('state.%s', nmstr);
            else
                statestr = 'state';
            end
            funcstr = func2str(fn);
            funcstr = regexp(funcstr, "\.([^.]*)$", 'tokens');
            funcstr = funcstr{1}{1};
            funcstr = horzcat(nmstr, {funcstr});
            funcstr = strjoin(funcstr, '.');
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

            if isempty(nmstr)
                state0str = 'state0';
            else
                state0str = sprintf('state0.%s', nmstr);
            end

            callstr = sprintf('%s = %s(%s, %s, dt);', statestr, funcstr, statestr, state0str);
            
        end

        function callstr = drivingForceFuncCallSetupFn(propfunction)

            [funcstr, statestr, nmstr] = propfunction.setupCallStringElements();
            
            callstr = sprintf('%s = %s(%s, drivingForces);', statestr, funcstr, statestr);
            
        end

        function callstr = literalFunctionCallSetupFn(propfunction)

        % First, we define string callstr1, which when evalulated will define the function called fn.
        %
        % Second, with the second string callstr2, the function fn is called with the correct argument (model and state
        % with the correct submodel setup).
            
            [~, statestr, nmstr] = propfunction.setupCallStringElements();

            callstr1 = sprintf('fn = %s;\n', func2str(propfunction.fn));

            if isempty(nmstr)
                nmstr = 'model';
            else
                nmstr = sprintf('model.%s', nmstr);
            end
            callstr2 = sprintf('%s = fn(%s, %s);', statestr, nmstr, statestr);

            callstr = [callstr1, callstr2];
            
        end


    end
    
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
