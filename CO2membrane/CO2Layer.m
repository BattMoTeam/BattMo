classdef CO2Layer < BaseModel

    properties

        Boundary
        
        constants
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

    end

    function model = registerVarAndPropfuncNames(model)

            varnames = {};
            %  fluxes [mol/s]
            varnames{end + 1} = VarName({}, 'Fluxes', nGas);
            %  pressures [pascal]
            varnames{end + 1} = VarName({}, 'Pressures', nGas);
            %  total pressure [pascal]
            varnames{end + 1} = 'Pressure';
            %  boundary sources [mol/s]
            varnames{end + 1} = VarName({}, 'BcSources', nGas);

            
            % mass conservation equations
            varnames{end + 1} = VarName({}, 'massConses', nGas);
            
            % Boundary condition flux definition equation
            varnames{end + 1} = VarName({'Boundary'}, 'bcDefinitions', nGas);

            % Boundary condition flux definition equation
            varnames{end + 1} = VarName({'Boundary'}, 'feedBcDefinitions', nGas);
            % Boundary condition flux definition equation
            varnames{end + 1} = VarName({'Boundary'}, 'permeateBcDefinitions', nGas);
            
            % Feed Control 
            varnames{end + 1} = VarName({'Boundary'}, 'feedControls', nGas);            
            % Permeate Control 
            varnames{end + 1} = VarName({'Boundary'}, 'permeateControls', nGas);
            
            model = model.registerVarNames(varnames);

    end
    
end
