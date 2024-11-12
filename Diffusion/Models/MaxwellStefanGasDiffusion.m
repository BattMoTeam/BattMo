classdef MaxwellStefanGasDiffusion < MaxwellStefanDiffusion
% Model for ideal gas

    properties

        referencePressure = 1*atm % use in chemical potential definition but should not influence the result as it
                                  % cancels when we take the gradient
        
    end
    
    methods

        function model = MaxwellStefanGasDiffusion(inputparams)

            model = model@MaxwellStefanDiffusion(inputparams);

        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@MaxwellStefanDiffusion(model);

            ncomp = model.numberOfComponents;

            fn = @MaxwellStefanGasDiffusion.updateDensity;
            inputvarnames = {'molarWeight', ...
                             'pressure'   , ...
                             'temperature'};
            model = model.registerPropFunction({'density', fn, inputvarnames});

            fn = @MaxwellStefanGasDiffusion.updateDensity;
            inputvarnames = {{'Boundary', 'molarWeight'}, ...
                             {'Boundary', 'pressure'}   , ...
                             {'Boundary', 'temperature'}};
            outputvarname = VarName({'Boundary'}, 'density');
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            for icomp = 1 : ncomp
                
                fn = @MaxwellStefanGasDiffusion.updateChemicalPotentials;
                outputvarname = VarName({}, 'chemicalPotentials', ncomp, icomp);
                % warning('check dependency')
                inputvarnames  = {'temperature', ...
                                  'density'    , ...
                                  'molarWeight', ...
                                  VarName({}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanGasDiffusion.updateChemicalPotentials;
                outputvarname = VarName({'Boundary'}, 'chemicalPotentials', ncomp, icomp);
                % warning('check dependency')
                inputvarnames  = {VarName({'Boundary'}, 'temperature'), ...
                                  VarName({'Boundary'}, 'density')    , ...
                                  VarName({'Boundary'}, 'molarWeight'), ...
                                  VarName({'Boundary'}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                                
            end
            
        end

            
        function state = updateDensity(model, state)

            state.density = model.computeDensity(state.pressure     , ...
                                                 state.molarWeight, ...
                                                 state.temperature);
            
            state.Boundary.density = model.computeDensity(state.Boundary.pressure     , ...
                                                          state.Boundary.molarWeight, ...
                                                          state.Boundary.temperature);
            
        end

        function density = computeDensity(model, pressure, molarWeight, temperature)
        % Ideal gas law

            R = PhysicalConstants.R;
            
            density = pressure./(R*temperature).*molarWeight;
            
        end

        function state = updateChemicalPotentials(model, state)

            state.chemicalPotentials          = model.computeChemicalPotential(state.pressure);
            state.Boundary.chemicalPotentials = model.computeChemicalPotential(state.Boundary.pressure);
            
        end

        function chemicalPotentials = computeChemicalPotential(model, density, massFractions, molarWeight, temperature)


            R = PhysicalConstants.R;
            
            p0    = model.referencePressure;
            mws   = model.molarWeights;
            ncomp = model.numberOfComponents;

            for icomp = 1 : ncomp

                pcomp = (massFractions{icomp}/mws{icomp}).*density.*(R.*temperature);
                
                chemicalPotentials{icomp} = R*temperature.*log(pcomp./p0);
                
            end
            
        end
        
    end
    
end

