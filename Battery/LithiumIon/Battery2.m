classdef Battery2 < PhysicalModel

    properties
        
        % Components
        Electrolyte
        NegativeElectrode
        PositiveElectrode
    end
    
    methods
        
        function model = Battery2(params)
        % Instantiate the components : Electrolyte, NegativeElectrode, PositiveElectrode
            
            model.G = params.G;
            
            model.Electrolyte       = params.Electrolyte;
            model.NegativeElectrode = params.NegativeElectrode;
            model.PositiveElectrode = params.PositiveElectrode;
            
        end
        
        
    end
    
end
