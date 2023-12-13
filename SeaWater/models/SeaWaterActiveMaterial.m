classdef SeaWaterActiveMaterial < BaseModel

    properties
        
        con = PhysicalConstants();
        chargeCarrierName
        stochElectron
        etaMax
    
    end
    
    methods
        
        function model = SeaWaterActiveMaterial(inputparams)

            model = model@BaseModel();
            fdnames = {'chargeCarrierName', ...
                       'stochElectron', ...
                       'etaMax'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end

        
        function model = registerVarAndPropfuncNames(model)
            
            varnames = {};
            % Temperature
            varnames{end + 1} = 'T'; 
            % Reaction rate (in mol/s/m^3)
            varnames{end + 1} = 'R'; 
            % electric potential from electrode 
            varnames{end + 1} = 'phiElectrode'; 
            % electric potential from electrolyte 
            varnames{end + 1} = 'phiElectrolyte'; 
            % concentration from electrolyte
            varnames{end + 1} = 'cElectrolyte'; 
            % value for eta (overpotential)
            varnames{end + 1} = 'eta'; 
            % value of ENernst coefficient
            varnames{end + 1} = 'ENernst'; 

            model = model.registerVarNames(varnames);

            fn = @() SeaWaterActiveMaterial.updateEta;
            inputnames = {'phiElectrolyte', 'phiElectrode', 'ENernst'};
            model = model.registerPropFunction({'eta', fn, inputnames});
            
        end
        
        
        function state = updateEta(model, state)          

            phiElectrolyte = state.phiElectrolyte;
            phiElectrode   = state.phiElectrode;
            ENernst        = state.ENernst;

            eta = phiElectrode - phiElectrolyte - ENernst;
            
            state.eta = eta;
        end            
        

    end
        
end
    
