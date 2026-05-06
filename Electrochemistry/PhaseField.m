classdef PhaseField < BaseModel
    
    properties
        
        
        %% Input parameters

        % Standard parameters

        N             % discretization parameter
        mobilityFunc  % mobility function
        dMobilityFunc % derivative of mobility function
        energyFunc    % Free energy function
        
        
        %% Helper structures
        
        A   % Collocation mapping matrices
        dA  % Collocation mapping matrices, first derivative
        ddA % Collocation mapping matrices, second derivative
        
    end
    
    methods
        
        function model = PhaseField(inputparams)
        %
        % ``inputparams`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %
            model = model@BaseModel();
            
            % fdnames = {'N'            , ... 
            %            'mobilityFunc' , ...
            %            'dMobilityFunc', ...
            %            'energyFunc'}

            % model = dispatchParams(model, inputparams, fdnames);

            % model = model.setupSpectralMethod();
            
        end

        function model = setupSpectralModel();

            % assign collocation mappings matrices (A, dA, ddA)
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Spectral decomposition coefficents for c
            varnames{end + 1} = 'coefC';
            % Spectral decomposition coefficents for w
            varnames{end + 1} = 'coefW';
            % collocation values for c
            varnames{end + 1} = 'c';
            % collocation values for w
            varnames{end + 1} = 'w';
            % collocation values for dC
            varnames{end + 1} = 'dC';
            % collocation values for dW
            varnames{end + 1} = 'dW';
            % collocation values for ddW
            varnames{end + 1} = 'ddC';
            % collocation values for ddW
            varnames{end + 1} = 'ddW';
            % equation for c
            varnames{end + 1} = 'eqC';
            % equation for w
            varnames{end + 1} = 'eqW';

            model = model.registerVarNames(varnames);

            fn = @PhaseField.updateC;
            model = model.registerPropFunction({'c', fn, {'coefC'}});
            
            fn = @PhaseField.updateW;
            model = model.registerPropFunction({'w', fn, {'coefW'}});

            fn = @PhaseField.updatedC;
            model = model.registerPropFunction({'dC', fn, {'coefC'}});
            
            fn = @PhaseField.updatedW;
            model = model.registerPropFunction({'dW', fn, {'coefW'}});

            fn = @PhaseField.updateddC;
            model = model.registerPropFunction({'ddC', fn, {'coefC'}});
            
            fn = @PhaseField.updateddW;
            model = model.registerPropFunction({'ddW', fn, {'coefW'}});
            
            fn = @PhaseField.updateEqC;
            model = model.registerPropFunction({'eqC', fn, {'c', 'dC', 'dW', 'ddW'}});

            fn = @PhaseField.updateEqW;
            model = model.registerPropFunction({'eqW', fn, {'w', 'ddC', 'c'}});
            
        end

        function model = setupForSimulation(model)
            
            model = model.equipModelForComputation();
            % model = model.setupScalings([]);
            % 
        end

        function state = updateC(model, state)

            state.c = model.A*state.coefC;
            
        end

        function state = updateW(model, state)
        end
        
        function state = updatedC(model, state)
        end
        
        function state = updatedW(model, state)
        end
        
        function state = updateddC(model, state)
        end
        
        function state = updateddW(model, state)
        end
        
        function state = updateEqC(model, state)
        end
        
        function state = updateEqW(model, state)
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
