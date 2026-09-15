classdef Function

% Parent class for an input function. The function can be given using different format. Each format is implemented with
% its own class. Use the separates function setupFunction to initiate the function

    properties
        
        functionFormat % String that can take values
                       % - 'tabulated',
                       % - 'string expression',
                       % - 'named function',
                       % - 'constant'

        argumentList % cell array of string which describes the expected argument of the function. Mainly for documentation at the moment or to infer number of argument

        %% helper

        numberOfArguments
        
    end

    methods


        function fn = Function(jsonstruct)

            fdnames = {'functionFormat', ...
                       'argumentList'};

            fn = dispatchParams(fn, jsonstruct, fdnames);

            fn.numberOfArguments = numel(fn.argumentList);

        end

        function y = eval(fn, varargin)
            
            error('virtual function');
            
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
