classdef SwellingCoating < Coating
% Implementation of model presented in :
% @article{Chandrasekaran_2010,
%   author =       {Chandrasekaran, Rajeswari and Magasinski, Alexandre and Yushin, Gleb and Fuller, Thomas F.},
%   title =        {Analysis of Lithium Insertion/Deinsertion in a Silicon Electrode Particle at Room Temperature},
%   year =         2010,
%   journal =      {Journal of The Electrochemical Society},
%   url =          {http://dx.doi.org/10.1149/1.3474225},
%   publisher =    {The Electrochemical Society},
% }
%
% For this model, we assume for the moment no binder and conductive additive.
%
    properties
        
        molarMass
        referenceFillInLevel % fill-in value which corresponds to the given volumeFraction of the coating

        includeHydrostaticStress = false;

        %% Computed at initialization
        
        zeroFillInVolumeFraction
        maximumTotalConcentration % The concentration is taken over the total volume (('total concentration') =
                                  % ('volume fraction')*('concentration in active material'))
        
    end
    
    methods
        
        function model = SwellingCoating(inputparams)

            model = model@Coating(inputparams)
            fdnames = {'molarMass', ...
                       'referenceFillInLevel'};
            
            model = dispatchParams(model, inputparams, fdnames);

            %% we assume for the moment no binder and conductive additive.
            % 
            if ~(all(model.volumeFractions == [1; 0; 0]))
                error('For the swelling model, the current implementation requires that there is no binder and conductive additive. To remove those, set their mass fractions equal to zero.');
            end
            
            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            
            %% Compute zero fill-in volume fraction from reference fill-in and given volume fraction
            %
            
            compmodel = model;
            compmodel = compmodel.registerVarAndPropfuncNames();
            compmodel = compmodel.removePropFunction({am, sd, 'x'});
            compmodel = compmodel.setupComputationalGraph();

            clear state
            state.(am).(sd).x = model.referenceFillInLevel;
            state = compmodel.evalVarName(state, {am, sd, 'radiusElongation'});

            re = state.(am).(sd).radiusElongation;

            model.zeroFillInVolumeFraction = model.volumeFraction/(re^3);

            %% Dispatch value to submodel
            model.(am).(sd).zeroFillInVolumeFraction = model.zeroFillInVolumeFraction;

            %% Compute maximum Lithium concentration with respect to total volume
            
            compmodel.zeroFillInVolumeFraction = model.zeroFillInVolumeFraction;
            
            clear state
            state.(am).(sd).x = 1;
            state = compmodel.evalVarName(state, 'volumeFraction');
            vfmax = state.volumeFraction;
            
            model.maximumTotalConcentration = vfmax*model.(am).(sd).saturationConcentration;

        end

        %% Declaration of the Dynamical Variables and Function of the model (setup of varnameList and propertyFunctionList)

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@Coating(model);
            
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {'volumeFraction'              , ...
                        {am, itf, 'volumetricSurfaceArea'}};

            model = model.registerVarNames(varnames);

            if model.includeHydrostaticStress
                model = model.registerVarName('hydrostaticStress');
            end
            
            fn = @SwellingCoating.updateVolumetricSurfaceArea;
            model = model.registerPropFunction({{am, itf, 'volumetricSurfaceArea'}, fn, {{am, sd, 'radiusElongation'}}});
            
            fn =  @SwellingCoating.updateRvol;
            model = model.registerPropFunction({{am, sd, 'Rvol'}, fn, {{am, itf, 'intercalationFlux'}, {am, itf, 'volumetricSurfaceArea'}}});

            fn =  @SwellingCoating.updateSolidIntercalationFlux;
            model = model.registerPropFunction({{am, sd, 'intercalationFlux'}, fn, {{am, itf, 'intercalationFlux'}}});

            if model.includeHydrostaticStress
                
                fn = @SwellingCoating.updateHydrostaticStress;
                model = model.registerPropFunction({'hydrostaticStress', fn, {{am, sd, 'cAverage'}, {am, itf, 'cElectrodeSurface'}}});

                fn  = @SwellingCoating.updateReactionRate;
                inputnames = {{am, itf, 'T'}, {am, itf, 'j0'}, {am, itf, 'eta'}, {am, sd, 'radius'}, 'hydrostaticStress'};
                model = model.registerPropFunction({{am, itf, 'intercalationFlux'}, fn, inputnames});
                
            end
            
            fn  = @SwellingCoating.updateConductivity;
            model = model.registerPropFunction({'conductivity', fn, {'volumeFraction'}});

            fn = @ActiveMaterial.updateVolumeFraction;
            model = model.registerPropFunction({'volumeFraction', fn, {{am, sd, 'radiusElongation'}}});
            
        end

        function newstate = addVariablesAfterConvergence(model, newstate, state)

            newstate = addVariablesAfterConvergence@Coating(model, newstate, state);
            newstate.volumeFraction = state.volumeFraction;

        end
        
        function state = updateVolumetricSurfaceArea(model, state)
            
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';


            delta = state.(am).(sd).radiusElongation;
            
            state.(am).(itf).volumetricSurfaceArea =  model.(am).(itf).volumetricSurfaceArea./delta;

        end

        function state = updateRvol(model, state)

            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            vsa = state.(am).(itf).volumetricSurfaceArea;
            R   = state.(am).(itf).intercalationFlux;
            
            Rvol = vsa.*R;

            state.(am).(sd).Rvol = Rvol;
            
        end


        function state = updateSolidIntercalationFlux(model, state)
            
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            state.(am).(sd).intercalationFlux = state.(am).(itf).intercalationFlux;
            
        end
        
        function state = updateConductivity(model, state)

            brugg = model.bruggemanCoefficient;

            vf = state.volumeFraction;
            
            % setup effective electrical conductivity using Bruggeman approximation
            state.conductivity = model.electronicConductivity.*vf.^brugg;

        end

        function state = updateReactionRate(model, state)
        % Same as in the interface class but uses the Butler Volmer
        % equation including stress

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            n     = model.(am).(itf).numberOfElectronsTransferred;
            F     = model.(am).(itf).constants.F;
            alpha = model.(am).(itf).chargeTransferCoefficient;

            T     = state.(am).(itf).T;
            j0    = state.(am).(itf).j0;
            eta   = state.(am).(itf).eta;
            sigma = state.hydrostaticStress;
            
            R = ButlerVolmerEquation_withStress(j0, alpha, n, eta, sigma, T);

            state.(am).(itf).intercalationFlux = R/(n*F); % reaction rate in mol/(s*m^2)

        end


        function state = updateCurrent(model, state)
            
            sigma = state.conductivity;
            phi   = state.phi;

            j = assembleFlux(model, phi, sigma);

            state.j = j;
            
        end
        
        function state = updateVolumeFraction(model, state)

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            vf = model.zeroFillInVolumeFraction;
            
            delta = state.(am).(sd).radiusElongation;
            
            state.volumeFraction = vf*delta.^3;
            
        end

        
        function state = updateHydrostaticStress(model, state)

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';

            E = 0;
            %Uncomment the expression of E  above for taking into account the stress
            %E         = 1e+11;
            nu    = 0.27;
            Omega = 4.25e-06;

            cSurface  = state.(am).(itf).cElectrodeSurface;
            cAverage  = state.(am).(sd).cAverage;

            sigma = ( (2.* E.* Omega)/(9.*(1-nu)) ) .* (cAverage - cSurface);

            state.hydrostaticStress = sigma;
            
        end

        %% Useful Functions
        function molarVolumeLithiated = computeMolarVolumeLithiated(model, theta)
        % cf equation 2 in [ref1]

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';

            molarVolumeSi = 1.2e-05;
            molarVolumeLi = 9e-06;

            molarVolumeLithiated = (4/15)*(molarVolumeSi + 3.75*theta*molarVolumeLi);
            
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
