classdef SingleCellElectrolyte < BaseModel
    
    methods
        
        function model = SingleCellElectrolyte(paramobj)
            model = model@BaseModel();
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Electrolyte potential
            varnames{end + 1} = 'phi';
            % Mass conservation equation (sum of ion flux  from cathode and anode vanishes)
            varnames{end + 1} = 'massCons';

            model = model.registerVarNames(varnames);
            
        end
        
    end
    
end
