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
        
        

        function model = ActiveMaterial(inputparams)
        %
        % ``inputparams`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %
            model = model@BaseModel();
            
            fdnames = {'N'            , ... 
                       'mobilityFunc' , ...
                       'dMobilityFunc', ...
                       'energyFunc'}

            model = dispatchParams(model, inputparams, fdnames);

            
        end
        
        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {'coefC'};
            varnames = {'coefW'};
            varnames = {'c'};
            varnames = {'w'};
            varnames = {'dC'};
            varnames = {'dW'};
            varnames = {'ddC'};
            varnames = {'ddW'};
            
            model = model.registerVarNames(varnames);

                
            % fn = @ActiveMaterial.updateControl;
            % fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            % model = model.registerPropFunction({'I', fn, {}});
                

        end

        function model = setupForSimulation(model)
            
            model = model.equipModelForComputation();
            % model = model.setupScalings([]);
            % 
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
