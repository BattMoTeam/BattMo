classdef Battery2 < PhysicalModel

    properties
        
        % Components
        Electrolyte
        NegativeElectrode
        PositiveElectrode
    
        coupTerms
    end
    
    methods
        
        function model = Battery2(params)
        % Instantiate the components : Electrolyte, NegativeElectrode, PositiveElectrode
            
            model.G = params.G;
            
            model.Electrolyte       = params.Electrolyte;
            model.NegativeElectrode = params.NegativeElectrode;
            model.PositiveElectrode = params.PositiveElectrode;
            
        end
    
        function state = setupElectrodeElectrolyteCoupling(model, state)
        % shortcuts:
        % elyte : Electrolyte
        % ne : NegativeElectrode
        % pe : PositiveElectrode
            
        end
        
        
        function state = updateActiveMaterialFromElyte(model, state)
        % shortcuts:
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
        % eac   : ElectrodeActiveComponent
        % am    : ActiveMaterial
            
            
        end
        
        
    end
    
end
