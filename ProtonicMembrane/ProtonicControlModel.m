classdef ProtonicControlModel < ProtonicMembraneGasSupplyBc
    
    methods
        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@ProtonicMembraneGasSupplyBc(model);

            nGas = model.nGas;
            
            varnames = {};
            
            % Rate (volumetric)
            % NOTE : Rate > 0 when injecting (opposity convention as massFlux)
            varnames{end + 1} = 'rate';           
            % Pressure equation for the control pressure
            varnames{end + 1} = 'pressureEq';
            % Rate equation for the control equation
            varnames{end + 1} = VarName({}, 'massFluxEqs', nGas);
            % control equation 
            varnames{end + 1} = 'setupEq';

            model = model.registerVarNames(varnames);

            model = model.removeVarName(VarName({}, 'boundaryEquations', nGas));

            for igas = 1 : nGas
                fn = @ProtonicControlModel.updateMassFractions;
                inputnames = {VarName({}, 'volumefractions', nGas, igas), ...
                              'density', ...
                              'rate'};
                outputname = VarName({}, 'massFluxes', nGas, igas);
                model = model.registerPropFunction({outputname, fn, inputnames});
            end
            
        end

        function state = updateMassFractions(model, state)

            nGas = model.nGas;
            
            vfs  = state.volumefractions;
            rho  = state.density;
            rate = state.rate;

            for igas = 1 : nGas

                % Note the minus sign below due to different sign convention for rate and massFluxes, respectively
                % inwards and outwards.
                state.massFluxes{igas} = - vfs{igas}.*rho.*rate;
                
            end
            
        end
        
    end
    
end
