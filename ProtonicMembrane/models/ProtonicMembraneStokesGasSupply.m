classdef ProtonicMembraneStokesGasSupply < ProtonicMembraneGasSupply
    
    properties


    end
    
    methods
        
        function model = ProtonicMembraneStokesGasSupply(inputparams)
            
            model = model@ProtonicMembraneGasSupply(inputparams);

        end


        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@ProtonicMembraneGasSupply(model);
            
            nGas = model.nGas;

            varnames = {};
            % Velocities / m/s (the 2D spatial indexing follows from the stokes solver setup)
            varnames = {'veloticies'};
            % Stokes equation
            varnames = {'stokesEquation'};
            % Total fluxes at the grid faces
            varnames = {'totalFlux'};
            
            model = model.registerVarNames(varnames);

            fn = @ProtonicMembraneStokesGasSupply.setupStokesEquation;
            inputvarnames = {'velocities', 'pressure'};
            model = model.registerPropFunction({'stokesEquation', fn, inputvarnames});
            
            fn = @ProtonicMembraneStokesGasSupply.updateTotalFlux;
            inputvarnames = {'velocities', 'pressure'};
            model = model.registerPropFunction({'totalFlux', fn, inputvarnames});
            
            for igas = 1 : nGas

                fn = @ProtonicMembraneStokesGasSupply.setupMassFluxes;
                inputvarnames = {VarName({}, 'massfractions', nGas, igas), ...
                                 VarName({}, 'densities', nGas, igas), ...
                                 'density', 'totalFlux'};
                outputvarname = VarName({}, 'massFluxes', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            end
            
        end


        
    end
    
end

