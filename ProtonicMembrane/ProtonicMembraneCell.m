classdef ProtonicMembraneCell < BaseModel
    
    properties
        
        % Temperature
        T
        % Structure with physical constants
        constants

        Anode
        Cathode
        Electrolyte
        
    end
    
    methods
        
        function model = ProtonicMembraneCell(paramobj)

            model = model@BaseModel();

            % model.operators = localSetupOperators(model.G);

            model.constants = PhysicalConstants();

            model.Anode       = ProtonicMembraneElectrode(paramobj.Anode);
            model.Cathode     = ProtonicMembraneElectrode(paramobj.Cathode);
            model.Electrolyte = ProtonicMembraneElectrolyte(paramobj.Electrolyte);
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % Total Hydrogen current
            varnames{end + 1} = 'I';
            % Cell potential
            varnames{end + 1} = 'E';

            model = model.registerVarNames(varnames);
        
            fn = @ProtonicMembraneCell.dispatchE;
            inputnames = {{elyte, 'E'}};
            model = model.registerPropFunction({{an, 'E'}, fn, inputnames});
            model = model.registerPropFunction({{ct, 'E'}, fn, inputnames});
            
            fn = @ProtonicMembraneCell.setupHpSources;
            inputnames = {{an, 'iBV'}, {ct, 'iBV'}};
            model = model.registerPropFunction({{elyte, 'sourceHp'}, fn, inputnames});
            
            fn = @ProtonicMembraneCell.setupElSources;
            inputnames = {{an, 'E'}, {ct, 'E'}, {elyte, 'E'}};
            model = model.registerPropFunction({{elyte, 'sourceEl'}, fn, inputnames});
                        
        end
        
        
    end
    
end
