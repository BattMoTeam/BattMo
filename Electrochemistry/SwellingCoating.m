classdef SwellingCoating < Coating

    % Same class as ActiveMaterial but with a new primaryVariable
    % (porosity). Some properties depending on porosity are thus no more constant

    
    properties
        
        molarMass
        
    end
    
    methods
        
        function model = SwellingCoating(inputparams)

            model = model@Coating(inputparams)
            fdnames = {'molarMass'};
            
            model = dispatchParams(model, inputparams, fdnames);

        end

        %% Declaration of the Dynamical Variables and Function of the model (setup of varnameList and propertyFunctionList)

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@Coating(model);
            
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {'porosity'                    , ...
                        'volumeFraction'              , ...
                        'porosityAccum'               , ...
                        'porositySource'              , ...
                        'porosityFlux'                , ...
                        'volumeCons'                  , ...
                        'hydrostaticStress'           , ...
                        {am, itf, 'volumetricSurfaceArea'}};

            model = model.registerVarNames(varnames);
            
            fn = @SwellingCoating.updatePorosityAccum;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'porosityAccum', fn, {'porosity','volumeFraction'}});

            fn = @SwellingCoating.updatePorositySource;
            inputparams = {{am, itf, 'intercalationFlux'}, {am, itf, 'volumetricSurfaceArea'}, {am, sd, 'cAverage'}, {am, sd,'radius'}};
            model = model.registerPropFunction({'porositySource', fn, inputparams});

            fn = @SwellingCoating.updatePorosityFlux;
            model = model.registerPropFunction({'porosityFlux', fn, {}});
            
            fn = @SwellingCoating.updateVolumeConservation;
            model = model.registerPropFunction({'volumeCons', fn, {'porosityAccum', 'porositySource', 'porosityFlux'}});
            
            fn = @SwellingCoating.updateHydrostaticStress;
            model = model.registerPropFunction({'hydrostaticStress', fn, {{am, sd, 'cAverage'}, {am, itf, 'cElectrodeSurface'}}});
            
            fn = @SwellingCoating.updateVolumeFraction;
            model = model.registerPropFunction({{'volumeFraction'}, fn, {'porosity'}});
            model = model.registerPropFunction({{am, sd, 'volumeFraction'}, fn, {'porosity'}});

            fn = @SwellingCoating.updateVolumetricSurfaceArea;
            model = model.registerPropFunction({{am, itf, 'volumetricSurfaceArea'}, fn, {{am, sd, 'radius'}, {am, sd, 'volumeFraction'}}});
            fn =  @SwellingCoating.updateRvol;
            model = model.registerPropFunction({{am, sd, 'Rvol'}, fn, {{am, itf, 'intercalationFlux'}, {am, itf, 'volumetricSurfaceArea'}}});

            fn  = @SwellingCoating.updateReactionRate;
            inputnames = {{am, itf, 'T'}, {am, itf, 'j0'}, {am, itf, 'eta'}, {am, sd, 'radius'}, 'hydrostaticStress'};
            model = model.registerPropFunction({{am, itf, 'intercalationFlux'}, fn, inputnames});
            
            fn  = @SwellingCoating.updateReactionRateCoefficient;
            inputnames = {{am, itf, 'cElectrolyte'}, {am, itf, 'cElectrodeSurface'}, {am, sd, 'radius'}, {am, itf, 'T'}};
            model = model.registerPropFunction({{am, itf, 'j0'}, fn, inputnames});

            fn  = @SwellingCoating.updateVolumeFraction;
            model = model.registerPropFunction({'volumeFraction', fn, {'porosity'}});
            model = model.registerPropFunction({{am, sd, 'volumeFraction'}, fn, {'porosity'}});

            fn  = @SwellingCoating.updateConductivity;
            model = model.registerPropFunction({'conductivity', fn, {'volumeFraction'}});
            
        end
        

        function state = updateRvol(model, state)

            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            vsa = state.(am).(itf).volumetricSurfaceArea;
            R   = state.(am).(itf).R;
            
            Rvol = vsa.*R;

            state.(am).(sd).Rvol = Rvol;
            
        end

        function state = updateConductivity(model, state)

            brugg = model.bruggemanCoefficient;

            vf = state.volumeFraction;
            
            % setup effective electrical conductivity using Bruggeman approximation
            state.conductivity = model.electronicConductivity.*vf.^brugg;

        end

        function state = updateCurrent(model, state)
        % Assemble electrical current which is stored in :code:`state.j`

            sigma = state.conductivity;
            phi   = state.phi;

            j = assembleFlux(model, phi, sigma);

            state.j = j;

        end

        function state = updateReactionRateCoefficient(model, state)

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            if model.(am).(itf).useJ0Func

                error('not checked');
                
                computeJ0 = model.(am).(itf).computeJ0Func;
                cmax      = model.(am).(itf).cmax;
                
                c = state.(am).(itf).cElectrodeSurface;
                R = state.(am).(sd).radius;

                theta = c/cmax;
                
                j0 = computeJ0(theta);
                
            else
                
                Tref = 298.15;  % [K]

                cmax          = model.(am).(itf).saturationConcentration;
                k0            = model.(am).(itf).reactionRateConstant;
                Eak           = model.(am).(itf).activationEnergyOfReaction;
                n             = model.(am).(itf).numberOfElectronsTransferred;
                F             = model.(am).(itf).constants.F;
                R             = model.(am).(itf).constants.R;
                R_delithiated = model.(am).(sd).particleRadius;

                T      = state.(am).(itf).T;
                cElyte = state.(am).(itf).cElectrolyte;
                c      = state.(am).(itf).cElectrodeSurface;
                radius = state.(am).(sd).radius;
                
                % Calculate reaction rate constant
                k = k0.*exp(-Eak./R.*(1./T - 1/Tref));

                %k = k.*(R_delithiated./radius).^2;

                % We use regularizedSqrt to regularize the square root function and avoid the blow-up of derivative at zero.
                th = 1e-3*cmax;
                coef = cElyte.*(cmax - c).*c;
                coef(coef < 0) = 0;
                
                j0 = k.*regularizedSqrt(coef, th).*n.*F;
                
            end
            
            state.(am).(itf).j0 = j0;

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

            r  = state.(am).(sd).radius;
            r0 = model.(am).(sd).particleRadius; %.* (r0./r).^2

            state.(am).(itf).R = R/(n*F); % reaction rate in mol/(s*m^2)

        end

        
        function state = updateVolumeFraction(model, state)

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            amf = model.volumeFractions(model.compInds.(am));
            
            porosity = state.porosity;

            vf = 1 - porosity;

            state.volumeFraction           = vf;
            state.(am).(sd).volumeFraction = amf*vf;
            
        end


        function state = updateVolumetricSurfaceArea(model, state)
        % Geometric result giving the volumetric surface area (cf bottom of
        % page 3 in Modelling capacity fade in silicon-graphite composite electrodes for
        % lithium-ion batteries Shweta Dhillon, Guiomar Hernández, Nils P. Wagner, Ann Mari Svensson,
        % Daniel Brandell ([ref1])

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';

            vf     = state.(am).(sd).volumeFraction;
            radius = state.(am).(sd).radius;

            state.(am).(itf).volumetricSurfaceArea = (3.*vf)./radius;
            
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

        %% Implementation of a new equation : the volume conservation Equation
        % Reference : eq 2 in [ref1]
        
        function state = updatePorosityAccum(model, state, state0, dt)

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';

            state.porosityAccum = (state.porosity - state0.porosity)./dt;
            
        end
        
        function state = updatePorositySource(model, state)
        % cf eq 2 in [ref1]
            
            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';

            cmax   = model.(am).(itf).saturationConcentration;
            theta0 = model.(am).(itf).guestStoichiometry0;
            
            R  = state.(am).(itf).R;
            a  = state.(am).(itf).volumetricSurfaceArea;       
            c  = state.(am).(sd).cAverage;

            theta = c/cmax;
            
            molarVolumeLithiated   = model.computeMolarVolumeLithiated(theta);
            molarVolumeDelithiated = model.computeMolarVolumeLithiated(theta0);
            
            state.porositySource = 0 .* a.*R.*(molarVolumeLithiated - molarVolumeDelithiated);
            
        end

        function state = updatePorosityFlux(model, state)
        % No Flux term (hack to create one)

            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            D = 0.*state.porosity;
            c = 0.*state.porosity;
            porosityFlux = assembleFlux(model, c, D);
            
            state.porosityFlux = porosityFlux;
            
        end

        function state = updateVolumeConservation(model, state)
            
            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            itf = 'Interface';

            flux   = state.porosityFlux;
            source = state.porositySource;
            accum  = state.porosityAccum;

            cons = assembleConservationEquation(model, flux, 0, source, accum);
            
            state.volumeCons = cons;
            
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
