classdef CO2captureChannelBoundary < BaseModel
    
    properties

        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        
    end
    
    methods
        
        function model = CO2captureChannelBoundary(gasSpecies)
            
            model = CO2capture.setupGasStructures(model, gasSpecies);
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            nGas = model.nGas;
            
            varnames = {};
            %  fluxes [mol/s] (positive sign for outer fluxes)
            varnames{end + 1} = VarName({}, 'rates', nGas);
            %  mol fractions
            varnames{end + 1} = VarName({}, 'molFractions', nGas);
            %  total pressure [pascal]
            varnames{end + 1} = 'pressure';
            %  total flux
            varnames{end + 1} = 'totalRate';
            %  total rate equation
            varnames{end + 1} = 'totalRateEquation';
            %  mol fraction equation (sum of mol fraction should be equal to one)
            varnames{end + 1} = 'molFractionConstraintEq';
            %  mol fraction equations. Provides the mol fraction at the boundary either from control value or from upwinding
            varnames{end + 1} = VarName({}, 'molFractionEquations', nGas);
            %  control equation, sets either the prescribed rate or pressure values
            varnames{end + 1} = 'controlEquation';
            
            model = model.registerVarNames(varnames);

            fn = @CO2captureChannelBoundary.totalRateEquation;
            inputvarnames = {VarName({}, 'rates', nGas), ...
                             'totalRate'};
            model = model.registerPropFunction({'totalRateEquation', fn, inputvarnames});

            fn = @CO2captureChannelBoundary.totalRateEquation;
            inputvarnames = {VarName({}, 'molFractions', nGas)};
            model = model.registerPropFunction({'molFractionConstraintEq', fn, inputvarnames});            
            
        end


    end
    
end

