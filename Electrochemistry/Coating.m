classdef Coating < ElectronicComponent

    properties

        %% Sub-Models

        ActiveMaterial
        Binder
        ConductingAdditive

        %% Input Parameters

        % Standard parameters
        density              % the mass density of the material (symbol: rho)
        bruggemanCoefficient % the Bruggeman coefficient for effective transport in porous media (symbol: beta)
        
        % Advanced parameters (used if given, otherwise computed)
        effectiveThermalConductivity
        volumeFraction
        
        %% Computed parameters at model setup

        activematerial_case % simple_particle or double_particle
        massFractions       % mass fractions of each components
        compInds            % index of the sub models in the massFractions structure
        
    end

    methods
        
        function model = Coating(paramobj)

            model = model@ElectronicComponent(paramobj);
            
            fdnames = {'density'                     , ...
                       'bruggemanCoefficient'        , ...
                       'effectiveThermalConductivity', ...
                       'volumeFraction'};

            model = dispatchParams(model, paramobj, fdnames);

            model.ActiveMaterial     = ActiveMaterial(paramobj.ActiveMaterial);
            model.Binder             = Binder(paramobj.Binder);
            model.ConductingAdditive = ConductingAdditive(paramobj.ConductingAdditive);
            
        end
        
    end
end

