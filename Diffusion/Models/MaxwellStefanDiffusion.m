classdef MaxwellStefanDiffusion < BaseModel
% Maxwell stefan diffusion model as provided by
% @article{curtiss1999multicomponent,
%   title={Multicomponent diffusion},
%   author={Curtiss, CF and Bird, R Byron},
%   journal={Industrial \& Engineering Chemistry Research},
%   volume={38},
%   number={7},
%   pages={2515--2522},
%   year={1999},
%   publisher={ACS Publications}
% }
%
% We consider (for now) isothermal case and no external forces
    
    properties

        % Diffusion maxtrix coefficients of dimension NxN (where N is the number of components)
        diffusionMatrix
        % Cell array with the component names
        compNames
        % Struct which maps a component name (provided as a struct field name) to component index (value returned by the
        % struct)
        compInds
        % Number of components
        numberOfComponents
        % Molecular weights (array with N values, where N is the number of components)
        molecularWeights

        % Boundary sub-model
        Boundary

        externalCouplingTerms % Cell array with external coupling terms

        % helpers
        bcfaces % Boundary faces that are included. For the boundary faces that are not included we have no flux boundary conditions
        
    end
    
    methods

        function model = MaxwellStefanDiffusion(inputparams)

            model = model@BaseModel();
            
            fdnames = {'G'               , ...
                       'diffusionMatrix' , ...
                       'compNames'       , ...
                       'molecularWeights', ...
                       'externalCouplingTerms'};
            
            model = dispatchParams(model, inputparams, fdnames);

            model.numberOfComponents = numel(model.compNames);

            model.subModelNameList = {'Boundary'};
            
            ncomp = model.numberOfComponents;
            
            for icomp = 1 : ncomp

                name = model.compNames{icomp};
                compInds.(name) = icomp;
                
            end

            model.compInds = compInds;
            
            % Check symmetry and the properties 2.8 in main reference

        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            ncomp = model.numberOfComponents;

            % Define first the variables that are common to the main domain and the boundary
            varnames = {};
            % Chemical potential
            varnames{end + 1} = VarName({}, 'chemicalPotentials', ncomp);
            % Diffusional driving forces
            varnames{end + 1} = VarName({}, 'diffusionForces', ncomp);            
            % Diffusion fluxes
            varnames{end + 1} = VarName({}, 'diffusionFluxes', ncomp);
            % Mass fractions
            varnames{end + 1} = VarName({}, 'massFractions', ncomp);
            % Mass fraction equation. The mass fractions sum to one.
            varnames{end + 1} = 'massFractionEquation';
            % Pressure
            varnames{end + 1} = 'pressure';
            % Temperature
            varnames{end + 1} = 'temperature';
            % Total density
            varnames{end + 1} = 'density';            
            % Concentration
            varnames{end + 1} = 'concentration';
            
            model = model.registerVarNames(varnames);
            
            for ivar = 1 : numel(varnames)
                varname = varnames{ivar};
                if ~isa(varname, 'VarName')
                    varname = VarName({'Boundary'}, varname);
                else
                    varname.namespace = {'Boundary'};
                end
                varnames{ivar} = varname;
            end
            
            model = model.registerVarNames(varnames);

            % Register the variables specific the boundary

            varnames = {};
            % Definition equations for the diffusion forces at the boundary
            varnames{end + 1} = VarName({'Boundary'}, 'diffusionForceEquations', ncomp);
            model = model.registerVarNames(varnames);
            
            % Register the variable specific to main domain (mass conservation equations)
            varnames = {};
            % Mass accumulation terms
            varnames{end + 1} = VarName({}, 'massAccums', ncomp);
            % Mass source terms
            varnames{end + 1} = VarName({}, 'sources', ncomp);
            % Boundary mass source terms
            varnames{end + 1} = VarName({}, 'boundarySources', ncomp);
            % Mass conservation equations
            varnames{end + 1} = VarName({}, 'massConses', ncomp);
            
            model = model.registerVarNames(varnames);
            
            model = model.setAsStaticVarName('temperature');
            model = model.setAsStaticVarName({'Boundary', 'temperature'});
            
            fn = @MaxwellStefanDiffusion.updateConcentration;
            inputvarnames = {'density', VarName({}, 'massFractions', ncomp)};
            outputvarname = VarName({}, 'concentration');
            model = model.registerPropFunction({'concentration', fn, inputvarnames});
            
            inputvarnames = {VarName({'Boundary'}, 'density'), ...
                             VarName({'Boundary'}, 'massFractions', ncomp)};
            outputvarname = VarName({'Boundary'}, 'concentration');
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
            fn = @MaxwellStefanDiffusion.updateDiffusionFluxes;
            inputvarnames = {'density'                          , ...
                             VarName({}, 'massFractions', ncomp), ...
                             VarName({}, 'diffusionForces', ncomp)};
            outputvarnames = VarName({}, 'diffusionFluxes', ncomp);
            model = model.registerPropFunction({outputvarnames, fn, inputvarnames});

            fn = @MaxwellStefanDiffusion.updateDiffusionFluxes;
            inputvarnames = {VarName({'Boundary'}, 'density')             , ...
                             VarName({'Boundary'}, 'massFractions', ncomp), ...
                             VarName({'Boundary'}, 'diffusionForces', ncomp)};
            outputvarnames = VarName({'Boundary'}, 'diffusionFluxes', ncomp);
            model = model.registerPropFunction({outputvarnames, fn, inputvarnames});
            
            fn = @MaxwellStefanDiffusion.updateMassFractionEquation;
            inputvarnames = {VarName({}, 'massFractions', ncomp)};
            model = model.registerPropFunction({'massFractionEquation', fn, inputvarnames});

            fn = @MaxwellStefanDiffusion.updateMassFractionEquation;
            inputvarnames = {VarName({'Boundary'}, 'massFractions', ncomp)};
            outputvarname = VarName({'Boundary'}, 'massFractionEquation');
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                        
            for icomp = 1 : ncomp
                
                fn = @MaxwellStefanDiffusion.updateDiffusionForces;
                outputvarname = VarName({}, 'diffusionForces', ncomp, icomp);
                inputvarnames  = {'temperature'                                  , ...
                                  'pressure'                                     , ...
                                  'density'                                      , ...
                                  'concentration'                                , ...
                                  VarName({}, 'chemicalPotentials', ncomp, icomp), ...
                                  VarName({}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanDiffusion.updateBoundaryDiffusionForcesEquations;
                outputvarname = VarName({'Boundary'}, 'diffusionForceEquations', ncomp, icomp);
                inputvarnames  = {'temperature'                                            , ...
                                  'pressure'                                               , ...
                                  'density'                                                , ...
                                  'concentration'                                          , ...
                                  VarName({}, 'chemicalPotentials', ncomp, icomp)          , ...
                                  VarName({}, 'massFractions', ncomp, icomp)               , ...
                                  VarName({'Boundary'}, 'temperature')                     , ...
                                  VarName({'Boundary'}, 'pressure')                        , ...
                                  VarName({'Boundary'}, 'density')                         , ...
                                  VarName({'Boundary'}, 'concentration')                   , ...
                                  VarName({'Boundary'}, 'chemicalPotentials', ncomp, icomp), ...
                                  VarName({'Boundary'}, 'massFractions', ncomp, icomp)     , ...
                                  VarName({'Boundary'}, 'diffusionForces', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
                fn = @MaxwellStefanDiffusion.updateMassAccums;
                fn = {fn, @(propfunc) PropFunction.accumFuncCallSetupFn(propfunc)};
                outputvarname = VarName({}, 'massAccums', ncomp, icomp);
                inputvarnames  = { 'density', ...
                                   VarName({}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanDiffusion.updateBoundarySources;
                outputvarname = VarName({}, 'boundarySources', ncomp, icomp);
                inputvarnames  = {VarName({'Boundary'}, 'diffusionFluxes', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanDiffusion.updateMassConses;
                outputvarname = VarName({}, 'massConses', ncomp, icomp);
                inputvarnames  = {VarName({}, 'massAccums'     , ncomp, icomp), ...
                                  VarName({}, 'diffusionFluxes', ncomp, icomp), ...
                                  VarName({}, 'sources'        , ncomp, icomp), ...
                                  VarName({}, 'boundarySources', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});                
                
            end
            
        end
        
    end
    
end

