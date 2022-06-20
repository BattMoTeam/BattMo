classdef ProtonicMembraneElectrode < BaseModel
    
    properties
        
        T  % Temperature
        constants
        
        % coefficient in Buttler-Volmer
        beta
        % Exchange current density
        iBV_0
        % Limiting current densities
        anodeLimitingCurrentDensity
        cathodeLimitingCurrentDensity
        % charge
        z
        
    end
    
    methods
        
        function model = ProtonicMembraneElectrode(paramobj)

            model = model@BaseModel();

            fdnames = {'G', ...
                       'beta', ...
                       'iBV_0', ...
                       'anodeLimitingCurrentDensity', ...
                       'cathodeLimitingCurrentDensity', ...
                       'z'};
            
            % model.operators = localSetupOperators(model.G);

            model.constants = PhysicalConstants();
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % Reaction rate for Buttler-Volmer
            varnames{end + 1} = 'iBV';
            % Open circuit potential
            varnames{end + 1} = 'Eocp';
            % Electrolyte electronic potential
            varnames{end + 1} = 'E';
            % Electrolyte electronic potential
            varnames{end + 1} = 'phiElectrolyte';
            
            model = model.registerVarNames(varnames);
        
            fn = @ProtonicMembraneElectrode.updateButtlerVolmerRate;
            inputnames = {'Eocp', 'E'};
            model = model.registerPropFunction({'iBV', fn, inputnames});
        
        end
        
        
        function state = updateButtlerVolmerRate(model, state)
            
            R  = model.constants.R;
            F  = model.constants.F;
            T  = model.T;
            z  = model.z;
            ia = model.anodeLimitingCurrentDensity;
            ic = model.cathodeLimitingCurrentDensity;
            i0 = model.iBV_0;
            
            Eocp = state.Eocp
            E    = state.E;
            
            f = F*z/(RT);
            eta = E - Eocp
            
            iBV = i0*(exp(-beta*f*eta) - exp((1 - beta)*f*eta))/(1 + (i0/ic)*exp(-beta*f*eta) - (i0/ia)*exp(-(1 - beta)*f*eta));
            
            state.iBV = iBV;
            
        end
        
        
        
        
        
    end
    
end
