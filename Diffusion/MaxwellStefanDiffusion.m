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
        
    end
    
    methods

        function model = MaxwellStefanDiffusion(inputparams)

            model = model@BaseModel();
            
            fdnames = {'G'                 , ...
                       'diffusionMatrix'   , ...
                       'compNames'         , ...
                       'numberOfComponents', ...
                       'molecularWeights'};
            
            model = dispatchParams(model, inputparams, fdnames);

            model.numberOfComponents = numel(model.compNames);

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

            varnames = {};
            % Chemical potential
            varnames{end + 1} = VarName({}, 'chemicalPotentials', ncomp);
            % Diffusional driving forces
            varnames{end + 1} = VarName({}, 'diffusionForces', ncomp);            
            % Diffusion fluxes
            varnames{end + 1} = VarName({}, 'diffusionFluxes', ncomp);
            % Mass fractions
            varnames{end + 1} = VarName({}, 'massFractions', ncomp);
            % Pressure
            varnames{end + 1} = 'pressure';
            % Temperature
            varnames{end + 1} = 'temperature';
            % Total density
            varnames{end + 1} = 'density';            
            % concentration
            varnames{end + 1} = 'concentration';
            
            model = model.registerVarNames(varnames);

            model = model.setAsStaticVarName('temperature');
            
            fn = @MaxwellStefanDiffusion.updateConcentration;
            inputvarnames = {'density', VarName({}, 'massFractions', ncomp)};
            model = model.registerPropFunction({'concentration', fn, inputvarnames});

            fn = @MaxwellStefanDiffusion.updateDiffusionFluxes;
            inputvarnames = {'density'                          , ...
                             VarName({}, 'massFractions', ncomp), ...
                             VarName({}, 'diffusionForces', ncomp)};
            outputvarnames = VarName({}, 'diffusionFluxes', ncomp);
            model = model.registerPropFunction({outputvarnames, fn, inputvarnames});

            fn = @MaxwellStefanDiffusion.updateDiffusionForces;
            for icomp = 1 : ncomp
                outputvarname = VarName({}, 'diffusionForces', ncomp, icomp);
                inputvarnames  = {'temperature'                                  , ...
                                  'pressure'                                     , ...
                                  'density'                                      , ...
                                  'concentration'                                , ...
                                  VarName({}, 'chemicalPotentials', ncomp, icomp), ...
                                  VarName({}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            end
            
        end
        
    end
    
end

