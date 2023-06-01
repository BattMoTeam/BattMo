classdef SwellingMaterial < ActiveMaterial

    % Same class as ActiveMaterial but with a new primaryVariable :
    % porosity and many properties which are no more constants but
    % dependant on the porosity
    
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
                        'volumeCons'};

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
            rp    = state.radius;
            vsf   = state.Interface.volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
            
            eqs = {};
            eqs{end + 1} = state.chargeCons;
            eqs{end + 1} = scalingcoef*state.(sd).massCons;
            eqs{end + 1} = scalingcoef*state.(sd).solidDiffusionEq;
            eqs{end + 1} = state.volumeCons
            
            names = {'chargeCons', ...
                     'massCons', ...
                     'solidDiffusionEq', ...
                     'volumeCons'};
            
            types = {'cell', 'cell', 'cell', 'cell'};

            primaryVars = model.getPrimaryVariables();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end


%% Update of the variables which were already defined in ActiveMaterial class but depending on the volumeFraction or the volumetric
% surface area which are no more constant parameters

        %% Same as in Active Material but for a non constant volumeFraction    
        function state = assembleAccumTerm(model, state, state0, dt)
        % Used when diffusionModelType == 'simple'
            
            vols   = model.G.cells.volumes;
            amFrac = model.activeMaterialFraction;

            vf     = state.volumeFraction;
            c  = state.c;
            c0 = state0.c;

            state.massAccum = vols.*vf.*amFrac.*(c - c0)/dt;
            
        end

        %% Same as in Active Material but for a non constant volumetricSurfaceArea    
        function state = updateRvol(model, state)

            vsa = state.Interface.volumetricSurfaceArea;
            
            Rvol = vsa.*state.Interface.R;
            state.Rvol = Rvol;

            if model.use_particle_diffusion
                state.SolidDiffusion.Rvol = Rvol;
            end
            
        end
        
        %% Same as in Active Material but for a non constant volumeFraction    
        function state = updateAverageConcentration(model, state)

            sd  = 'SolidDiffusion';

            am_frac  = model.activeMaterialFraction;
            vols     = model.G.cells.volumes;

            vf       = state.volumeFraction;
            c = state.(sd).cAverage;

            vols = am_frac*vf.*vols;

            cAverage = sum(c.*vols)/sum(vols);

            state.cAverage = cAverage;
            
        end

        %% Same as in Active Material but for a non constant volumeFraction    
        function state = updateSOC(model, state)

            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            am_frac  = model.activeMaterialFraction;
            vols     = model.G.cells.volumes;
            cmax     = model.(itf).cmax;
            theta100 = model.(itf).theta100;
            theta0   = model.(itf).theta0;

            vf       = state.volumeFraction;
            c        = state.(sd).cAverage;

            theta = c/cmax;
            m     = (1 ./ (theta100 - theta0));
            b     = -m .* theta0;
            SOC   = theta*m + b;
            vol   = am_frac*vf.*vols;
            
            SOC = sum(SOC.*vol)/sum(vol);

            state.SOC = SOC;
            
        end
   
%% Update of the new variables (which are constant parameters in the case of ActiveMaterial)


        function state = updateRadius(model, state, state0)

            radius_0  = model.SolidDiffusion.rp;
            densitySi = model.density;
            cmaxLi    = model.Interface.cmax;
            
            c = state.SolidDiffusion.cAverage;

            
            MolarMassSi   = model.molarMass;
            molarVolumeSi = MolarMassSi/densitySi;
            molarVolumeLi = 8.8 * 1E-6;

            % We cannot anymore express <x> as a ratio of between the current and maximal concentrations of Li because the radius of the 
            % particle is changing. We have to express it as a ratio of matter quantities N/Nmax. The expression above is from a
            % rearrangement of equations 11 and 14 in 'Analysis of Lithium Insertion/Deinsertion in a Silicon Electrode
            % Particle at Room Temperature Rajeswari
            % Chandrasekaran,Alexandre Magasinski, Gleb Yushin, and Thomas F. Fuller'*


            % Only delithiation radius evolution is coded here. Schould be
            % easy to imlement for lithiation (cf equation 11 in the same
            % paper)

            c_ratio = c/cmaxLi;

            ratio_delith = (3.75.*molarVolumeLi)./(molarVolumeSi+3.75.*molarVolumeLi);

            radius = radius_0 .* ((1-ratio_delith) ./ ((1./3.8) + ratio_delith.*c_ratio)) .^ (1/3);

            state.radius = radius;
            
            if model.use_particle_diffusion
                state.SolidDiffusion.radius = radius;
            end
            
        end

        function state = updateVolumeFraction(model, state)

            porosity = state.porosity;
            vf = 1 - porosity;

            state.volumeFraction = vf;
            state.Interface.volumeFraction = vf;

             if model.use_particle_diffusion
                state.SolidDiffusion.volumeFraction = vf;
             end
             
        end
               
        function state = updateVolumetricSurfaceArea(model, state)
            
            amf = model.activeMaterialFraction;

            vf     = state.Interface.volumeFraction;
            radius = state.radius;

            vsa = (3.*vf.*amf)./radius;

            state.Interface.volumetricSurfaceArea = vsa;

            if model.use_particle_diffusion
                
                amf    = model.activeMaterialFraction;

                vf     = state.Interface.volumeFraction;
                radius = state.SolidDiffusion.radius;

                vsa = (3.*vf.*amf)./radius;

                state.SolidDiffusion.volumetricSurfaceArea = vsa;
                
            end
        end

        function state = updateConductivity(model, state)
            
            brugg = model.BruggemanCoefficient;

            porosity = state.porosity;

            vf = 1 - porosity;
                   
            % setup effective electrical conductivity using Bruggeman approximation
            state.conductivity = model.electricalConductivity.*vf.^brugg;
            
        end

        function state = updateReactionRateCoefficient(model, state)
            if model.Interface.useJ0Func

                computeJ0 = model.Interface.computeJ0Func;
                cmax      = model.Interface.cmax;
                theta0    = model.Interface.theta0;
                theta100  = model.Interface.theta100;
                radius0   = model.rp;
                
                c = state.Interface.cElectrodeSurface;
                radius = state.radius;

                %cmin = theta0*cmax;
                %cmax = theta100*cmax;

                % We cannot anymore express the state of charge as ratio of concentrations because the radius of the particle is
                % changing. We have to express it  in terms of quantities
                % of matter: x = N/Nmax

                scalingTerm = (radius./radius0) .^3;

                x = (c./cmax) .* scalingTerm;

                soc = (x - xmin)./(xmax - xmin);
                
                j0 = computeJ0(soc);

            else
                
                Tref = 298.15;  % [K]

                cmax = model.Interface.cmax;
                k0   = model.Interface.k0;
                Eak  = model.Interface.Eak;
                n    = model.Interface.n;
                F    = model.Interface.constants.F;
                R    = model.Interface.constants.R;
                radius_0 = model.SolidDiffusion.rp;

                T      = state.Interface.T;
                cElyte = state.Interface.cElectrolyte;
                c      = state.Interface.cElectrodeSurface;
                radius = state.SolidDiffusion.radius;


               % Necessary to define a c', cf eq 9 paper Analysis of Lithium Insertion/Deinsertion in a Silicon Electrode
               % Particle at Room Temperature Rajeswari Chandrasekaran,Alexandre Magasinski, Gleb Yushin,cand
               % Thomas F. Fuller
               
               scalingCoeff = (radius./radius_0).^3;
               c = c .* scalingCoeff;
 
                  
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




   %% Implementation of a new equation : the volume conservation Equation
   % ref : eq 2 in 'Modelling capacity fade in silicon-graphite composite electrodes for
   % lithium-ion batteries Shweta Dhillon, Guiomar Hernández, Nils P. Wagner, Ann 
   % Mari Svensson, Daniel Brandell
   
        function state = updatePorosityAccum(model, state, state0, dt)
            vols = model.G.cells.volumes;
            
            state.porosityAccum = vols.*(state.porosity - state0.porosity)./dt;
            
        end
            
        function state = updatePorositySource(model, state)
            vols = model.G.cells.volumes;
            
            molarVolumeLithiated = model.updateMolarVolumeLithiated(state);
            densitySi            = model.Interface.density;
            molarMassSi          = model.molarMass;

            a = state.Interface.volumetricSurfaceArea;       
            R = state.Interface.R;
            
            molarVolumeSi = molarMassSi/densitySi;

            state.porositySource = a.*R.*(molarVolumeLithiated - (4/15)*molarVolumeSi).*vols;
            
        end

         function state = updatePorosityFlux(model, state)
             %no Source. Definition schould be change to something cleaner
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

  


   %% Useful Function giving the molar volume of the lithiated silicon depending on its state of lithiation     
        function molarVolumeLitihated = updateMolarVolumeLithiated(model, state)
            
            c = state.SolidDiffusion.cAverage;

            densitySi = model.Interface.density;
            molarMassSi = model.molarMass;
            cmaxLi = model.Interface.cmax;

            molarVolumeSi = molarMassSi/densitySi;
            molarVolumeLi = 8.8 * 1E-6;
            

            molarVolumeLitihated = (4/15)*(molarVolumeSi + 3.75*(c/cmaxLi)*molarVolumeLi);
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
