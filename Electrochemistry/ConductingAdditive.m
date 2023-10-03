classdef ConductingAdditive < BaseModel

    properties

        %% Input Parameters
        
        %  Standard parameters
        
        electronicConductivity % the electronic conductivity of the material (symbol: sigma)
        density                % the mass density of the material (symbol: rho)
        massFraction           % the ratio of the mass of the material to the total mass of the phase or mixture (symbol: gamma)
        thermalConductivity    % Thermal conductivity of current collector
        specificHeatCapacity   % Heat capacity of current collector
        
    end

    methods
        
        function model = ConductingAdditive(paramobj)
            
            model = model@BaseModel();
            
            fdnames = {'electronicConductivity', ... 
                       'density'               , ...                
                       'massFraction'          , ...
                       'thermalConductivity'   , ...
                       'specificHeatCapacity'};

            model = dispatchParams(model, paramobj, fdnames);

        end
        
    end
end

