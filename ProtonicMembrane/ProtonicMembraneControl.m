classdef ProtonicMembraneControl < BaseModel
    
    properties
        
    end
    
    methods
        
        function model = ProtonicMembraneControl(paramobj)

            model = model@BaseModel();

        end

        function model = registerVarAndPropfuncNames(model)

            varnames = {};

            % Potential
            varnames{end + 1} = 'U';
            % Current
            varnames{end + 1} = 'I';            
            % Control Equation
            varnames{end + 1} = 'controlEquation';

            model = model.registerVarNames(varnames);
            
            fn = @ProtonicMembraneCell.setupHpSources;
            inputnames = {'U', 'I'};
            model = model.registerPropFunction({'controlEquation', fn, inputnames});
            
        end

        
    end
    
end
