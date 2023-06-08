classdef SwellingMaterial < ActiveMaterial

    % Same class as ActiveMaterial but with a new primaryVariable
    % (porosity). Some properties depending on porosity are thus no more constant

    
    properties

    end
    
    methods
        
        function model = SwellingMaterial(paramobj)
        % cf 'if isSwellingMaterial' in the constructor defintion of Active Material
            model = model@ActiveMaterial(paramobj)
        end

        %% Declaration of the Dynamical Variables and Function of the model (setup of varnameList and propertyFunctionList)
            % Implementation of all the new variables (primaryvariable
            % porosity and variables depending on porosity)
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ActiveMaterial(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {'porosity'       ,...
                        'radius'         ,...
                        'volumeFraction' ,...
                        'porosityAccum'  ,...
                        'porositySource' ,...
                        'porosityFlux'   ,...
                        'volumeCons'     ,...
                        'hydrostaticStress'};

            model = model.registerVarNames(varnames);
      

            fn = @SwellingMaterial.updatePorosityAccum;
            model = model.registerPropFunction({'porosityAccum', fn, {'porosity'}});

            fn = @SwellingMaterial.updatePorositySource;
            model = model.registerPropFunction({'porositySource', fn, {{itf, 'R'}, {itf, 'volumetricSurfaceArea'}}});

             fn = @SwellingMaterial.updatePorosityFlux;
            model = model.registerPropFunction({'porosityFlux', fn, {}});

            fn = @SwellingMaterial.updateVolumeConservation;
            model = model.registerPropFunction({'volumeCons', fn, {'porosityAccum', 'porositySource', 'porosityFlux'}});


            fn = @SwellingMaterial.updateConductivity;
            model = model.registerPropFunction({'conductivity', fn, {'volumeFraction'} });


            fn = @SwellingMaterial.updateHydrostaticStress
            model = model.registerPropFunction({'hydrostaticStress', fn, {sd,'cAverage'} });

            if model.use_particle_diffusion
                
                fn = @SwellingMaterial.updateRadius;
                model = model.registerPropFunction({'radius', fn, {{sd, 'cAverage'}}});
                model = model.registerPropFunction({{sd, 'radius'}, fn, {{sd, 'cAverage'}} });
                
                fn = @SwellingMaterial.updateVolumeFraction;
                model = model.registerPropFunction({{'volumeFraction'}, fn, {'porosity'}});
                model = model.registerPropFunction({{sd, 'volumeFraction'}, fn, {'porosity'}});
                model = model.registerPropFunction({{itf, 'volumeFraction'}, fn, {'porosity'}});


                fn = @SwellingMaterial.updateVolumetricSurfaceArea;
                model = model.registerPropFunction({{itf, 'volumetricSurfaceArea'}, fn, {{'radius'}, {itf, 'volumeFraction'}}});
                model = model.registerPropFunction({{sd, 'volumetricSurfaceArea'}, fn, {{sd,'radius'},{sd, 'volumeFraction'}}});

                fn =  @SwellingMaterial.updateRvol;
                model = model.registerPropFunction({{'Rvol'}, fn, {{itf,'R'}, {itf, 'volumetricSurfaceArea'}}});
                model = model.registerPropFunction({{sd, 'Rvol'}, fn, {{itf,'R'}, {itf, 'volumetricSurfaceArea'}}});

 

                
            else

                fn = @SwellingMaterial.updateRadius;
                model = model.registerPropFunction({'radius', fn, {{sd, 'cAverage'}}});
                
                fn = @ActiveMaterial.updateVolumeFraction;
                model = model.registerPropFunction({'volumeFraction', fn, {'porosity'}});
                model = model.registerPropFunction({{itf, 'volumeFraction'}, fn, {'porosity'}});

                fn = @SwellingMaterial.updateVolumetricSurfaceArea;
                model = model.registerPropFunction({{itf, 'volumetricSurfaceArea'}, fn, {{'radius'}, {itf, 'volumeFraction'}}});

                fn =  @SwellingMaterial.updateRvol;
                model = model.registerPropFunction({{'Rvol'}, fn, {{itf,'R'}, {itf, 'volumetricSurfaceArea'}}});

                
            end          
                      
                   
        end

        %% Similar to ActiveMaterial but adding the updates for the new variables
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            time = state0.time + dt;
            state = model.initStateAD(state);

            state = updateControl(model, state, drivingForces);


            %Specific to swelling Materials
            state                = model.updateVolumeFraction(state);
            state                = model.updateConductivity(state)    ;     
           
           
            state                = model.updateStandalonejBcSource(state);
            state                = model.updateCurrent(state);
            state.SolidDiffusion = model.SolidDiffusion.updateMassAccum(state.SolidDiffusion, state0.SolidDiffusion, dt);
            state                = model.dispatchTemperature(state);
            state.SolidDiffusion = model.SolidDiffusion.updateDiffusionCoefficient(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.updateFlux(state.SolidDiffusion);
            state                = model.updateConcentrations(state);
            
            %Specific to swelling Materials
            state                = model.updateRadius(state,state0);
            state                = model.updateVolumetricSurfaceArea(state);
            state                = model.updatePorosityAccum(state, state0, dt);
            state                = model.updatePorositySource(state);
            state                = model.updateVolumeConservation(state);
            

            state                = model.updatePhi(state);
            state.Interface      = model.Interface.updateReactionRateCoefficient(state.Interface);
            state.Interface      = model.Interface.updateOCP(state.Interface);
            state.Interface      = model.Interface.updateEta(state.Interface);
            state.Interface      = model.Interface.updateReactionRate(state.Interface);
            state                = model.updateRvol(state);
            state                = model.updateCurrentSource(state);
            state                = model.updateChargeConservation(state);
            state.SolidDiffusion = model.SolidDiffusion.updateMassSource(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.assembleSolidDiffusionEquation(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.updateMassConservation(state.SolidDiffusion);
            
            %% Setup equations and add some scaling
            n     = model.(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F     = model.(sd).constants.F;
            vol   = model.operators.pv;
            rp    = model.radius;
            vsf   = model.Interface.volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
            
            eqs = {};
            eqs{end + 1} = state.chargeCons;
            eqs{end + 1} = scalingcoef*state.(sd).massCons;
            eqs{end + 1} = scalingcoef*state.(sd).solidDiffusionEq;
            eqs{end + 1} = state.volumeCons;
            
            names = {'chargeCons', ...
                     'massCons', ...
                     'solidDiffusionEq', ...
                     'volumeCons'};
            
            types = {'cell', 'cell', 'cell', 'cell'};

            primaryVars = model.getPrimaryVariables();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function model = setupDependentProperties(model)           

            amFrac    = model.activeMaterialFraction;
            model.volumeFraction = 1 - model.porosity - amFrac ;
            vf = model.volumeFraction;
            brugg = model.BruggemanCoefficient;
            
            % setup effective electrical conductivity using Bruggeman approximation
            model.EffectiveElectricalConductivity = model.electricalConductivity.*vf.^brugg;

            
            if model.use_interparticle_diffusion
                
                interDiff = model.InterDiffusionCoefficient;
                
                
                model.EffectiveDiffusionCoefficient = interDiff.*(vf).^brugg;
                
            end

            if model.use_thermal
                % setup effective thermal conductivity
                model.EffectiveThermalConductivity = model.thermalConductivity.*vf.^brugg;
                model.EffectiveVolumetricHeatCapacity = model.specificHeatCapacity.*vf.*model.density;
            end
            
        end


%% Update for variables already defined in ActiveMaterial but depending on the volumeFraction or the volumetric
% surface area which are no more constant parameters

        % Same as in Active Material but for a non constant volumeFraction    
        function state = assembleAccumTerm(model, state, state0, dt)
        % Used when diffusionModelType == 'simple'
            
            vols   = model.G.cells.volumes;
            vf     = state.volumeFraction;
            c  = state.c;
            c0 = state0.c;

            state.massAccum = vols.*vf.*(c - c0)/dt;
            
        end

        % Same as in Active Material but for a non constant volumetricSurfaceArea    
        function state = updateRvol(model, state)

            vsa = state.Interface.volumetricSurfaceArea;
            R   = state.Interface.R;
            
            Rvol = vsa .* R;

            state.Rvol = Rvol;
            
            if model.use_particle_diffusion
                state.SolidDiffusion.Rvol = Rvol;
            end           
        end
        
        % Same as in Active Material but for a non constant volumeFraction    
        function state = updateAverageConcentration(model, state)

            sd  = 'SolidDiffusion';

            vols     = model.G.cells.volumes;

            vf       = state.volumeFraction;
            c        = state.(sd).cAverage;

            vols = vf.*vols;

            cAverage = sum(c.*vols)/sum(vols);

            state.cAverage = cAverage;
            
        end

        function state = updateReactionRateCoefficient(model, state)
            if model.Interface.useJ0Func

                computeJ0       = model.Interface.computeJ0Func;
                cmax            = model.Interface.cmax;
                theta0          = model.Interface.theta0;
                theta100        = model.Interface.theta100;
                R_delithiated   = model.SolidDiffusion.rp;
                
                c = state.Interface.cElectrodeSurface;
                R = state.radius;

                soc = model.computeSOC(c);
                
                j0 = computeJ0(soc);
            else
                
                Tref = 298.15;  % [K]

                cmax = model.Interface.cmax;
                k0   = model.Interface.k0;
                Eak  = model.Interface.Eak;
                n    = model.Interface.n;
                F    = model.Interface.constants.F;
                R    = model.Interface.constants.R;
                R_delithiated = model.SolidDiffusion.rp;

                T      = state.Interface.T;
                cElyte = state.Interface.cElectrolyte;
                c      = state.Interface.cElectrodeSurface;
                radius = state.SolidDiffusion.radius;

                molarVolumeSi   = 1.2e-05;
                molarVolumeLi   = model.constants.molarVolumeLi;
                Q = (3.75.*molarVolumeLi)./(molarVolumeSi);


               % Necessary to rescale c and cmax, cf eq 9 in [ref3]               
               scalingCoeff = (radius./R_delithiated).^3;
               c = c .* scalingCoeff;
               cmax = (1+Q) * cmax;
                  
               % Calculate reaction rate constant
               k = k0.*exp(-Eak./R.*(1./T - 1/Tref));
   
               %   We use regularizedSqrt to regularize the square root function and avoid the blow-up of derivative at zero.
                th = 1e-3* model.Interface.cmax;
                coef = cElyte.*(cmax - c).*c;
                coef(coef < 0) = 0;
                j0 = k.*regularizedSqrt(coef, th).*n.*F;

                
            end
            
            state.Interface.j0 = j0;

        end


        function state = updateConductivity(model, state)
            
            brugg = model.BruggemanCoefficient;

            vf    = state.Interface.volumeFraction;
                   
            % setup effective electrical conductivity using Bruggeman approximation
            state.conductivity = model.electricalConductivity.*vf.^brugg;
            
        end

        function state = updateCurrentSource(model, state)
            
            F    = model.Interface.constants.F;
            vols = model.G.cells.volumes;
            n    = model.Interface.n;

            Rvol = state.Rvol;
            
            r = state.SolidDiffusion.radius;
            r0 = model.SolidDiffusion.rp;
            state.eSource = - vols.*Rvol*n*F; % C/s
            
        end

         function state = updateReactionRate(model, state)
        % Same as in the interface class but uses the Butler Volmer
        % zquation including stress
            n     = model.Interface.n;
            F     = model.Interface.constants.F;
            alpha = model.Interface.alpha;

            T   = state.Interface.T;
            j0  = state.Interface.j0;
            eta = state.Interface.eta;
            sigma = state.hydrostaticStress;
            
            R = ButlerVolmerEquation_withStress(j0, alpha, n, eta, sigma,T);

            state.Interface.R = R/(n*F); % reaction rate in mol/(s*m^2)

        end

   
%% Update of the new variables (variables which are constant parameters in the case of ActiveMaterial)


        function state = updateRadius(model, state)
     %eq 4 in Modelling capacity fade in silicon-graphite composite electrodes for lithium-ion batteries
     %Shweta Dhillon, Guiomar Hernández, Nils P. Wagner, Ann Mari Svensson, Daniel Brandell ([ref 1])

            R_delith = model.SolidDiffusion.rp;
            cmax     = model.Interface.cmax;

            cAverage = state.SolidDiffusion.cAverage;

            R_lith   = computeRadius(cAverage, cmax, R_delith);

            if R_lith < R_delith
                error('Radius inferior to R_delith, not coherent')
            end

            state.radius = R_lith;
            
            if model.use_particle_diffusion
                state.SolidDiffusion.radius = R_lith;
            end            
        end

        
        function state = updateVolumeFraction(model, state)

            porosity = state.porosity;
            amf = model.activeMaterialFraction;

            vf = 1 - porosity - amf;

            state.volumeFraction = vf;
            state.Interface.volumeFraction = vf;

             if model.use_particle_diffusion
                state.SolidDiffusion.volumeFraction = vf;
             end
             
        end


        function state = updateVolumetricSurfaceArea(model, state)
     % Geometric result giving the volumetric surface area (cf bottom of
     % page 3 in Modelling capacity fade in silicon-graphite composite electrodes for
     % lithium-ion batteries Shweta Dhillon, Guiomar Hernández, Nils P. Wagner, Ann Mari Svensson,
     % Daniel Brandell ([ref1])
            
            vf     = state.Interface.volumeFraction;
            radius = state.radius;

            vsa    = (3.*vf)./radius;

            state.Interface.volumetricSurfaceArea = vsa;

            if model.use_particle_diffusion

                vf     = state.Interface.volumeFraction;
                radius = state.SolidDiffusion.radius;

                vsa    = (3.*vf)./radius;

                state.SolidDiffusion.volumetricSurfaceArea = vsa;
                
            end
        end
      

        function state = updateHydrostaticStress(model, state)
            E = 0;
            %Uncomment the expression of E  above for taking into account the stress
            %E         = 1e+11;
            nu        = 0.27;
            Omega     = 4.25e-06;

            cSurface  = state.Interface.cElectrodeSurface;
            cAverage  = state.SolidDiffusion.cAverage;

            sigma = ( (2.* E.* Omega)/(9.*(1-nu)) ) .* (cAverage - cSurface);

            state.hydrostaticStress = sigma;
        end

       




   %% Implementation of a new equation : the volume conservation Equation
   % Reference : eq 2 in [ref1]
   
        function state = updatePorosityAccum(model, state, state0, dt)
            vols = model.G.cells.volumes;

            vf     = state.volumeFraction;
            
            state.porosityAccum =  vf .* vols.*(state.porosity - state0.porosity)./dt;
            
        end
            
        function state = updatePorositySource(model, state)
        % cf eq 2 in [ref1]
        
            vols   = model.G.cells.volumes;
            cmax = model.Interface.cmax;

            vf     = state.volumeFraction;
            c      = state.SolidDiffusion.cAverage;
            a      = state.Interface.volumetricSurfaceArea;       
            R      = state.Interface.R;

            molarVolumeLithiated   = model.updateMolarVolumeLithiated(c);
            molarVolumeDelithiated = model.updateMolarVolumeLithiated(0);

            r0 = model.SolidDiffusion.rp;
            r = computeRadius(c,cmax,r0);

            state.porositySource = a.*R.*(molarVolumeLithiated - molarVolumeDelithiated).*vols .* vf;
            
        end

         function state = updatePorosityFlux(model, state)
         %No source term
             D = 0.*state.volumeFraction;
             c = 0.*state.volumeFraction;
             porosityFlux = assembleFlux(model, c, D);
            
            state.porosityFlux = porosityFlux;
        end

        

        function state = updateVolumeConservation(model, state)
            
            flux = state.porosityFlux;
            source = state.porositySource;
            accum = state.porosityAccum;

            cons = assembleConservationEquation(model, flux, 0, source, accum);
            
            state.volumeCons = cons;
            
        end



  


   %% Useful Functions    
        function molarVolumeLitihated = updateMolarVolumeLithiated(model, c)
      % cf equation 2 in [ref1]

            molarVolumeSi   = 1.2e-05;
            molarVolumeLi   = model.constants.molarVolumeLi;

            soc = model.computeSOC(c);

            molarVolumeLitihated = (4/15)*(molarVolumeSi + 3.75*soc*molarVolumeLi);
        end

        
        function state = updateSOC(model, state)
        %not used

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


        function soc = computeSOC(model, c)
        % Useful fonction, cf expression of the soc in eq16 of [ref3]
            % shortcut
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


            m     = (1 ./ (theta100 - theta0));
            b     = -m .* theta0;
            soc   = theta*m + b;
            
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
