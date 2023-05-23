classdef SwellingMaterial < ActiveMaterial
    
    properties

    end
    
    methods
        
        function model = SwellingMaterial(paramobj)
        % ``paramobj`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %SameFunction as for ActiveMaterial but implemening
        %FullSolidDiffusionSwelling            
            model = model@ActiveMaterial(paramobj)
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ActiveMaterial(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {'radius'   , ...
                        'porosity' ,...
                        'volumeFraction',...
                        'porosityAccum',...
                        'porositySource',...
                        'porosityFlux' ,...
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


            fn = @SwellingMaterial.updateEffectiveElectricalConductivity;
            model = model.registerPropFunction({'EffectiveElectricalConductivity', fn, {'volumeFraction'} });


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

        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            time = state0.time + dt;
            state = model.initStateAD(state);

            state = updateControl(model, state, drivingForces);


            %new
            state                = model.updateVolumeFraction(state);
            state                = model.updateEffectiveElectricalConductivity(state)    ;     
           
           
            state                = model.updateStandalonejBcSource(state);
            state                = model.updateCurrent(state);
            state.SolidDiffusion = model.SolidDiffusion.updateMassAccum(state.SolidDiffusion, state0.SolidDiffusion, dt);
            state                = model.dispatchTemperature(state);
            state.SolidDiffusion = model.SolidDiffusion.updateDiffusionCoefficient(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.updateFlux(state.SolidDiffusion);
            state                = model.updateConcentrations(state);
            
            %new
            state                = model.updateRadius(state);
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

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@ActiveMaterial(model, state, problem, dx, drivingForces);
            
        end



        function state = updateConcentrations(model, state)

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            if model.use_particle_diffusion
            if strcmp(model.diffusionModelType, 'simple')
                state.(sd).cAverage = state.c;
            end
            
            state.(itf).cElectrodeSurface = state.(sd).cSurface;
            else
                state.(itf).cElectrodeSurface = state.c;
            end
            
        end

        function state = updateMassFlux(model, state)
        % Used when diffusionModelType == 'simple'

            D = model.EffectiveDiffusionCoefficient;
            
            c = state.c;

            massflux = assembleFlux(model, c, D);
            
            state.massFlux = massflux;

        end
            
        function state = assembleAccumTerm(model, state, state0, dt)
        % Used when diffusionModelType == 'simple'
            
            vols   = model.G.cells.volumes;
            vf     = state.volumeFraction;
            amFrac = model.activeMaterialFraction;

            c  = state.c;
            c0 = state0.c;

            state.massAccum = vols.*vf.*amFrac.*(c - c0)/dt;
            
        end

        function state = updateRvol(model, state)

            vsa = state.Interface.volumetricSurfaceArea;
            
            Rvol = vsa.*state.Interface.R;
            state.Rvol = Rvol;

            if model.use_particle_diffusion
                state.SolidDiffusion.Rvol = Rvol;
            end
            
        end
  
        function state = updateMassConservation(model, state)
        % Used when diffusionModelType == 'simple' or no particle diffusion
            
            flux = state.massFlux;
            source = state.massSource;
            accum = state.massAccum;
            
            cons = assembleConservationEquation(model, flux, 0, source, accum);
            
            state.massCons = cons;
            
        end

        function state = updateMassSource(model, state)
        % used when diffusionModelType == simple
            
            vols = model.G.cells.volumes;
            
            Rvol = state.Rvol;
            
            state.massSource = - Rvol.*vols;
            
        end
        
        function state = updateStandalonejBcSource(model, state)
            
            state.jBcSource = state.controlCurrentSource;

        end

        function state = updateCurrentSource(model, state)
            
            F    = model.Interface.constants.F;
            vols = model.G.cells.volumes;
            n    = model.Interface.n;

            Rvol = state.Rvol;
            
            state.eSource = - vols.*Rvol*n*F; % C/s
            
        end
        
        function state = updatePhi(model, state)
            state.Interface.phiElectrode = state.phi;
        end         
        
        function state = dispatchTemperature(model, state)
            state.Interface.T = state.T;
            state.SolidDiffusion.T = state.T;
        end

        function state = updatejBcSource(model, state)
            state.jBcSource = state.jCoupling + state.jExternal;
        end
        
        function state = updatejFaceBc(model, state)
            state.jFaceBc = state.jFaceCoupling + state.jFaceExternal;
        end
        
        function state = updatejExternal(model, state)
            state.jExternal = 0;
            state.jFaceExternal = 0;
        end

        function state = updatejCoupling(model, state)
            state.jCoupling = 0;
            state.jFaceCoupling = 0;
        end   
        
        function state = updateAverageConcentration(model, state)

            % shortcut
            sd  = 'SolidDiffusion';

            vf       = state.volumeFraction;
            am_frac  = model.activeMaterialFraction;
            vols     = model.G.cells.volumes;
            
            c = state.(sd).cAverage;

            vols = am_frac*vf.*vols;

            cAverage = sum(c.*vols)/sum(vols);

            state.cAverage = cAverage;
            
        end
        
        function state = updateSOC(model, state)

            % shortcut
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            vf       = state.volumeFraction;
            am_frac  = model.activeMaterialFraction;
            vols     = model.G.cells.volumes;
            cmax     = model.(itf).cmax;
            theta100 = model.(itf).theta100;
            theta0   = model.(itf).theta0;
            
            c = state.(sd).cAverage;

            theta = c/cmax;
            m     = (1 ./ (theta100 - theta0));
            b     = -m .* theta0;
            SOC   = theta*m + b;
            vol   = am_frac*vf.*vols;
            
            SOC = sum(SOC.*vol)/sum(vol);

            state.SOC = SOC;
            
        end
   

%Functions specific to swelling materials

        function state = updateRadius(model, state)

            radius_0  = model.SolidDiffusion.rp;
            densitySi = model.density;
            cmaxLi    = model.Interface.cmax;
            
            c = state.SolidDiffusion.cAverage;
            
            MolarMassSi   = 28.0855 * 1E-3;
            molarVolumeSi = MolarMassSi/densitySi;
            molarVolumeLi = 8.8 * 1E-6;
            
            radius = radius_0 * (1 + (3.75*molarVolumeLi*c)/(cmaxLi*molarVolumeSi))^(1/3);

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

        function state = updateEffectiveElectricalConductivity(model, state)
            
            brugg = model.BruggemanCoefficient;

            porosity = state.porosity;

            vf = 1 - porosity;
                   
            % setup effective electrical conductivity using Bruggeman approximation
            state.EffectiveElectricalConductivity = model.electricalConductivity.*vf.^brugg;
            
        end

        function state = updatePorosityAccum(model, state, state0, dt)
            
            state.porosityAccum = (state.porosity - state0.porosity)/dt;
            
        end
            
        function state = updatePorositySource(model, state)
            
            molarVolumeLithiated = model.updateMolarVolumeLithiated(state);
            densitySi            = model.Interface.density;

            a = state.Interface.volumetricSurfaceArea;       
            R = state.Interface.R;
            
            molarMassSi   = 28.0855 * 1E-3;
            molarVolumeSi = molarMassSi/densitySi;

            state.porositySource = -a.*R.*(molarVolumeLithiated - 3.75*molarVolumeSi);
            
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
                    
        function molarVolumeLitihated = updateMolarVolumeLithiated(model, state)

            sd  = 'SolidDiffusion';

            c = state.(sd).cAverage;

            densitySi = model.Interface.density;
            molarMassSi = 28.0855 * 1E-3;
            molarVolumeSi = molarMassSi/densitySi;

            molarVolumeLi = 8.8 * 1E-6;
            cmaxLi = model.Interface.cmax;

            molarVolumeLitihated = 4/15*(molarVolumeSi + 3.75*(c/cmaxLi)*molarVolumeLi);
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
