classdef FullSolidDiffusionSwellingModel < FullSolidDiffusionModel
% Implementation of model presented in :
% @article{Chandrasekaran_2010,
%   author =       {Chandrasekaran, Rajeswari and Magasinski, Alexandre and Yushin, Gleb and Fuller, Thomas F.},
%   title =        {Analysis of Lithium Insertion/Deinsertion in a Silicon Electrode Particle at Room Temperature},
%   year =         2010,
%   journal =      {Journal of The Electrochemical Society},
%   url =          {http://dx.doi.org/10.1149/1.3474225},
%   publisher =    {The Electrochemical Society},
% }

    properties

        referenceFillInLevel % Fill-in value which corresponds to the given radius

        %% Physical properties
        % Molar volumes. Values taken from reference, table I.
        
        molarVolumeSi = 1.2e-05; % m^3/mol
        molarVolumeLi = 9e-6; % m^3/mol
        
        %% Computed at initialization

        zeroFillInParticleRadius
        zeroFillInVolumetricSurfaceArea
        
        %% Advanced parameters
        % typically setup by swelling coating model
        zeroFillInVolumeFraction
        
    end

    methods

        function model = FullSolidDiffusionSwellingModel(inputparams)

            model = model@FullSolidDiffusionModel(inputparams);
            
            fdnames = {'referenceFillInLevel', ...
                       'zeroFillInVolumeFraction'};

            model = dispatchParams(model, inputparams, fdnames);

            %% compute zero fill in radius from reference fill-in and given particule radius

            compmodel = model;
            compmodel = compmodel.registerVarAndPropfuncNames();
            compmodel = compmodel.removePropFunction('x');
            compmodel = compmodel.setupComputationalGraph();

            clear state
            state.x = compmodel.referenceFillInLevel;
            state = compmodel.evalVarName(state, 'radiusElongation');

            re = state.radiusElongation;

            model.zeroFillInParticleRadius        = model.particleRadius/re;
            model.zeroFillInVolumetricSurfaceArea = model.volumetricSurfaceArea*re;
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@FullSolidDiffusionModel(model);

            varnames = {'radiusElongation', ...
                        'radius'          , ...
                        'intercalationFlux'};

            % Particle fill-in level (between 0 and 1). Depends only on stoichiometry
            % By definition, we have : x = (volumeFraction*concentration)/(volumeFraction_max*concentration_max)
            % Due to the special expression for the volume fraction, it reduces to a function of the stoichiometry only,
            % See paper in reference

            varnames{end + 1} = 'x';

            model = model.registerVarNames(varnames);

            model = model.unsetAsExtraVarName('cAverage');

            model = model.setAsExtraVarName('radius'); % the radius is not needed directly. We use the radius elongation value instead.
            
            fn = @FullSolidDiffusionSwellingModel.updateMassSource;
            model = model.registerPropFunction({'massSource', fn, {'intercalationFlux', 'radiusElongation'}});

            fn = @FullSolidDiffusionSwellingModel.updateRadius;
            model = model.registerPropFunction({'radius', fn, {'radiusElongation'}});

            fn = @FullSolidDiffusionSwellingModel.updateFillInLevel;
            model = model.registerPropFunction({'x', fn, {'cAverage'}});
            
            fn = @FullSolidDiffusionSwellingModel.updateRadiusElongation;
            model = model.registerPropFunction({'radiusElongation', fn, {'x'}});

            fn = @FullSolidDiffusionSwellingModel.updateMassAccum;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};            
            model = model.registerPropFunction({'massAccum', fn, {'c', 'radiusElongation'}});
            
        end

        function state = updateMassAccum(model, state, state0, dt)

            op = model.operators;

            c      = state.c;
            delta  = op.mapToParticle*state.radiusElongation;
            c0     = state0.c;
            delta0 = op.mapToParticle*state0.radiusElongation;

            state.massAccum = 1./(delta.*dt).*op.vols.*(delta.^3.*c - delta0.^3.*c0);
            
        end
        
        function state = updateRadius(model, state)

            state.radius = model.zeroFillInParticleRadius*state.radiusElongation;
            
        end
        
        function state = updateFillInLevel(model, state)

            cmax = model.saturationConcentration;
            Q    = (3.75.*model.molarVolumeLi)./(model.molarVolumeSi);

            cAverage = state.cAverage;

            theta = cAverage/cmax;

            state.x = theta./(1 + Q*(1 -theta));
            
        end
        
        function state = updateRadiusElongation(model, state)
        % eq 11 in ref paper

            x = state.x;

            Q = (3.75.*model.molarVolumeLi)./(model.molarVolumeSi);

            state.radiusElongation = (1 + Q.*x).^(1/3);
            
        end

        function state = updateMassSource(model, state)

            op  = model.operators;
            
            rp0  = model.zeroFillInParticleRadius;
            vf0  = model.zeroFillInVolumeFraction;
            vsa0 = model.zeroFillInVolumetricSurfaceArea;
            
            delta = state.radiusElongation;
            R     = state.intercalationFlux;

            massSource = - (R./(delta.^2)).*((4/3)*pi*rp0^3.*vsa0/vf0);
            massSource = op.mapFromBc*massSource;

            state.massSource = massSource;
            
        end

        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState.radiusElongation = state.radiusElongation;
            
        end 

    end
    
end


%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
    
