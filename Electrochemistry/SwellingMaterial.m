classdef SwellingMaterial < ActiveMaterial

    % Same class as ActiveMaterial but with a new primaryVariable
    % (porosity). Some properties depending on porosity are thus no more constant

    
    properties
        
        molarMass
        
    end
    
    methods
        
        function model = SwellingMaterial(paramobj)

            model = model@ActiveMaterial(paramobj)
            fdnames = {'molarMass'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
        end

        function model = setupDiffusionModel(model, paramobj)

            model.use_particle_diffusion = true;
            model.use_interparticle_diffusion = false;
            paramobj.SolidDiffusion.np = model.G.cells.num;
            model.SolidDiffusion = FullSolidDiffusionSwellingModel(paramobj.SolidDiffusion);
            
        end

        %% Declaration of the Dynamical Variables and Function of the model (setup of varnameList and propertyFunctionList)

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ActiveMaterial(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {'porosity'                    , ...
                        'volumeFraction'              , ...
                        'porosityAccum'               , ...
                        'porositySource'              , ...
                        'porosityFlux'                , ...
                        'volumeCons'                  , ...
                        'hydrostaticStress'           , ...
                        {itf, 'volumetricSurfaceArea'}};

            model = model.registerVarNames(varnames);
            
            fn = @SwellingMaterial.updatePorosityAccum;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'porosityAccum', fn, {'porosity'}});

            fn = @SwellingMaterial.updatePorositySource;
            model = model.registerPropFunction({'porositySource', fn, {{itf, 'R'}, {itf, 'volumetricSurfaceArea'}, {sd, 'cAverage'}}});

            fn = @SwellingMaterial.updatePorosityFlux;
            model = model.registerPropFunction({'porosityFlux', fn, {}});
            
            fn = @SwellingMaterial.updateVolumeConservation;
            model = model.registerPropFunction({'volumeCons', fn, {'porosityAccum', 'porositySource', 'porosityFlux'}});
            
            fn = @SwellingMaterial.updateConductivity;
            model = model.registerPropFunction({'conductivity', fn, {'volumeFraction'}});
                                               
            fn = @SwellingMaterial.updateHydrostaticStress;
            model = model.registerPropFunction({'hydrostaticStress', fn, {{sd, 'cAverage'}, {itf, 'cElectrodeSurface'}}});
            
            fn = @SwellingMaterial.updateVolumeFraction;
            model = model.registerPropFunction({{'volumeFraction'}, fn, {'porosity'}});
            model = model.registerPropFunction({{sd, 'volumeFraction'}, fn, {'porosity'}});

            fn = @SwellingMaterial.updateVolumetricSurfaceArea;
            model = model.registerPropFunction({{itf, 'volumetricSurfaceArea'}, fn, {{sd, 'radius'}, 'volumeFraction'}});

            fn =  @SwellingMaterial.updateRvol;
            model = model.registerPropFunction({'Rvol', fn, {{itf, 'R'}, {itf, 'volumetricSurfaceArea'}}});
            model = model.registerPropFunction({{sd, 'Rvol'}, fn, {{itf, 'R'}, {itf, 'volumetricSurfaceArea'}}});

            fn  = @SwellingMaterial.updateReactionRate;
            inputnames = {{itf, 'T'}, {itf, 'j0'}, {itf, 'eta'}, 'hydrostaticStress'};
            model = model.registerPropFunction({{itf, 'R'}, fn, inputnames});

            fn  = @SwellingMaterial.updateVolumeFraction;
            model = model.registerPropFunction({'volumeFraction', fn, {'porosity'}});
            model = model.registerPropFunction({{sd, 'volumeFraction'}, fn, {'porosity'}});
            
        end
        

        function state = updateRvol(model, state)
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            am_fraction  = model.activeMaterialFraction;

            vsa = state.(itf).volumetricSurfaceArea;
            R   = state.(itf).R;
            
            Rvol = vsa.*R;

            state.Rvol = Rvol;
            state.(sd).Rvol = Rvol;
            
        end
        
        % Same as in Active Material but for a non constant volumeFraction    
        function state = updateAverageConcentration(model, state)

            sd  = 'SolidDiffusion';

            vols     = model.G.cells.volumes;
            am_frac  = model.activeMaterialFraction;

            vf       = state.volumeFraction;
            c        = state.(sd).cAverage;

            vols = amf.*vf.*vols;

            cAverage = sum(c.*vols)/sum(vols);

            state.cAverage = cAverage;
            
        end

        function state = updateReactionRateCoefficient(model, state)
            
            if model.Interface.useJ0Func

                error('not checked');
                
                computeJ0       = model.Interface.computeJ0Func;
                cmax            = model.Interface.cmax;
                theta0          = model.Interface.theta0;
                theta100        = model.Interface.theta100;
                R_delithiated   = model.SolidDiffusion.rp;
                
                c = state.Interface.cElectrodeSurface;
                R = state.radius;

                soc = model.computeTheta(c);
                
                j0 = computeJ0(soc);
                
            else
                
                Tref = 298.15;  % [K]

                cmax          = model.Interface.cmax;
                k0            = model.Interface.k0;
                Eak           = model.Interface.Eak;
                n             = model.Interface.n;
                F             = model.Interface.constants.F;
                R             = model.Interface.constants.R;
                R_delithiated = model.SolidDiffusion.rp;

                T      = state.Interface.T;
                cElyte = state.Interface.cElectrolyte;
                c      = state.Interface.cElectrodeSurface;
                radius = state.SolidDiffusion.radius;

                molarVolumeSi = 1.2e-05;
                molarVolumeLi = model.constants.molarVolumeLi;
                
                Q = (3.75.*molarVolumeLi)./(molarVolumeSi);

                % Necessary to rescale c and cmax, cf eq 9 in [ref3]
                scalingCoeff = (radius./R_delithiated).^3;
                c            = c.*scalingCoeff;
                cmax         = (1 + Q)*cmax;
                
                % Calculate reaction rate constant
                k = k0.*exp(-Eak./R.*(1./T - 1/Tref));

                % k = k .* (R_delithiated./radius).^2;
                % We use regularizedSqrt to regularize the square root function and avoid the blow-up of derivative at zero.
                th = 1e-3* model.Interface.cmax;
                coef = cElyte.*(cmax - c).*c;
                coef(coef < 0) = 0;
                
                j0 = k.*regularizedSqrt(coef, th).*n.*F;
                
            end
            
            state.Interface.j0 = j0;

        end


        function state = updateConductivity(model, state)
            
            brugg = model.BruggemanCoefficient;

            vf    = state.volumeFraction;
            
            % setup effective electrical conductivity using Bruggeman approximation
            state.conductivity = model.electricalConductivity.*vf.^brugg;
            
        end

        function state = updateCurrentSource(model, state)
            
            F    = model.Interface.constants.F;
            vols = model.G.cells.volumes;
            n    = model.Interface.n;
            r0   = model.SolidDiffusion.rp;

            Rvol = state.Rvol;
            r    = state.SolidDiffusion.radius;

            % xavier : you do not use r?
            state.eSource = - vols.*Rvol*n*F; % C/s
            
        end

        function state = updateReactionRate(model, state)
        % Same as in the interface class but uses the Butler Volmer
        % equation including stress

            itf = 'Interface';
            
            n     = model.(itf).n;
            F     = model.(itf).constants.F;
            alpha = model.(itf).alpha;

            T     = state.(itf).T;
            j0    = state.(itf).j0;
            eta   = state.(itf).eta;
            sigma = state.hydrostaticStress;
            
            R = ButlerVolmerEquation_withStress(j0, alpha, n, eta, sigma, T);

            state.(itf).R = R/(n*F); % reaction rate in mol/(s*m^2)

        end

        
        function state = updateVolumeFraction(model, state)

            porosity = state.porosity;

            vf = 1 - porosity;

            state.volumeFraction                = vf;
            state.SolidDiffusion.volumeFraction = vf;
            
        end


        function state = updateVolumetricSurfaceArea(model, state)
        % Geometric result giving the volumetric surface area (cf bottom of
        % page 3 in Modelling capacity fade in silicon-graphite composite electrodes for
        % lithium-ion batteries Shweta Dhillon, Guiomar Hernández, Nils P. Wagner, Ann Mari Svensson,
        % Daniel Brandell ([ref1])

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            vf     = state.volumeFraction;
            radius = state.(sd).radius;

            vsa = (3.*vf)./radius;

            state.(itf).volumetricSurfaceArea = vsa;
            
        end
        

        function state = updateHydrostaticStress(model, state)

            E = 0;
            %Uncomment the expression of E  above for taking into account the stress
            %E         = 1e+11;
            nu    = 0.27;
            Omega = 4.25e-06;

            cSurface  = state.Interface.cElectrodeSurface;
            cAverage  = state.SolidDiffusion.cAverage;

            sigma = ( (2.* E.* Omega)/(9.*(1-nu)) ) .* (cAverage - cSurface);

            state.hydrostaticStress = sigma;
            
        end

        %% Implementation of a new equation : the volume conservation Equation
        % Reference : eq 2 in [ref1]
        
        function state = updatePorosityAccum(model, state, state0, dt)
            
            vols = model.G.cells.volumes;

            state.porosityAccum = vols.*(state.porosity - state0.porosity)./dt;
            
        end
        
        function state = updatePorositySource(model, state)
        % cf eq 2 in [ref1]
            
            vols   = model.G.cells.volumes;
            cmax   = model.Interface.cmax;
            theta0 = model.Interface.theta0;
            
            R  = state.Interface.R;
            a  = state.Interface.volumetricSurfaceArea;       
            c  = state.SolidDiffusion.cAverage;
            
            molarVolumeLithiated   = model.computeMolarVolumeLithiated(c);
            % Expression below schould be improved (cmin corresponding at theta0)
            molarVolumeDelithiated = model.computeMolarVolumeLithiated(1000);

            r0 = model.SolidDiffusion.rp;
            r  = computeRadius(c, cmax, r0);

            state.porositySource = a.*R.*(molarVolumeLithiated - molarVolumeDelithiated).*vols;
            
        end

        function state = updatePorosityFlux(model, state)
        % No Flux term (hack to create one)
            
            D = 0.*state.porosity;
            c = 0.*state.porosity;
            porosityFlux = assembleFlux(model, c, D);
            
            state.porosityFlux = porosityFlux;
            
        end

        function state = updateVolumeConservation(model, state)
            
            flux   = state.porosityFlux;
            source = state.porositySource;
            accum  = state.porosityAccum;

            cons = assembleConservationEquation(model, flux, 0, source, accum);
            
            state.volumeCons = cons;
            
        end


        %% Useful Functions    
        function molarVolumeLitihated = computeMolarVolumeLithiated(model, c)
        % cf equation 2 in [ref1]

            molarVolumeSi = 1.2e-05;
            molarVolumeLi = model.constants.molarVolumeLi;

            soc = model.computeTheta(c);

            molarVolumeLitihated = (4/15)*(molarVolumeSi + 3.75*soc*molarVolumeLi);
            
        end

        
        function state = updateSOC(model, state)

            error('should not be used')
            
            % shortcut
            itf = 'Interface';
            sd  = 'SolidDiffusion';


            vols     = model.G.cells.volumes;
            cmax     = model.(itf).cmax;
            theta100 = model.(itf).theta100;
            theta0   = model.(itf).theta0;
            R_delith = model.(sd).rp;
            molarVolumeSi = 1.2e-05;
            molarVolumeLi =  model.constants.molarVolumeLi;

            Q = (3.75.*molarVolumeLi)./(molarVolumeSi);
            
            c        = state.(sd).cAverage;
            vf       = state.volumeFraction;

            c_ratio = c/cmax;       

            R_lith = computeRadius(c,cmax,R_delith);

            theta = c_ratio .* ((R_lith./R_delith).^3) ./ (1+Q);


            m     = (1 ./ (theta100 - theta0));
            b     = -m .* theta0;
            SOC   = theta*m + b;
            vol   = vf.*vols;
            
            SOC = sum(SOC.*vol)/sum(vol);

            state.SOC = SOC;
            
        end


        function theta = computeTheta(model, c)
        % Useful fonction, cf expression of the soc in eq16 of [ref3]

            itf = 'Interface';
            sd  = 'SolidDiffusion';

            
            cmax     = model.(itf).cmax;
            theta100 = model.(itf).theta100;
            theta0   = model.(itf).theta0;
            R_delith = model.(sd).rp;
            molarVolumeSi = 1.2e-05;
            molarVolumeLi =  model.constants.molarVolumeLi;

            Q = (3.75.*molarVolumeLi)./(molarVolumeSi);
            
            c_ratio = c/cmax;       

            R_lith = computeRadius(c,cmax,R_delith);

            theta = c_ratio .* ((R_lith./R_delith).^3) ./ (1+Q);
            
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
