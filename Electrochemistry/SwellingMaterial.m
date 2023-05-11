classdef SwellingMaterial < ActiveMaterial
    
    properties

    end
    
    methods
        
        function model = SwellingMaterial(paramobj)
        % ``paramobj`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
            model = model@ActiveMaterial(paramobj);
            
            % hello
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ActiveMaterial(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {{sd, 'cAverage'}};
            model = model.registerVarNames(varnames);

            %Particle radius and thus porosity, volume fraction, volumetric
            %surface area and mass source are now variables.


            if model.use_particle_diffusion
                
                fn = @SwellingMaterial.updateRadius;
                model = model.registerPropFunction({'radius', fn, {{sd, 'cAverage'}}});
                model = model.registerPropFunction({{sd, 'radius'}, fn, {{sd, 'cAverage'}} });
                
                fn = @SwellingMaterial.updateVolumeFraction;
                model = model.registerPropFunction({{'volumeFraction'}, fn, {'porosity'}});
                model = model.registerPropFunction({{sd, 'volumeFraction'}, fn, {'porosity'}});
                model = model.registerPropFunction({{itf, 'volumeFraction'}, fn, {'porosity'}});


                fn = @SwellingMaterial.updateVolumetricSurfaceArea;
                model = model.registerPropFunction({{itf, 'volumetricSurfaceArea'}, fn, {'radius', {itf, 'volumeFraction'}}});
                model = model.registerPropFunction({{sd, 'volumetricSurfaceArea'}, fn, {{sd, 'radius'},{sd, 'volumeFraction'}}});

                
            else

                fn = @SwellingMaterial.updateRadius;
                model = model.registerPropFunction({'radius', fn, {{sd, 'cAverage'}}});
                
                fn = @ActiveMaterial.updateVolumeFraction;
                model = model.registerPropFunction({'volumeFraction', fn, {'porosity'}});
                model = model.registerPropFunction({{itf, 'volumeFraction'}, fn, {'porosity'}});

                fn = @SwellingMaterial.updateVolumetricSurfaceArea;
                model = model.registerPropFunction({{itf, 'volumetricSurfaceArea'}, fn, {{sd,'radius'}, {itf, 'volumeFraction'}}});

                
            end


            fn = @SwellingMaterial.updatePorosity;
            model = model.registerPropFunction({'porosity', fn, {}});

            fn = @SwellingMaterial.updateEffectiveElectricalConductivity;
            model = model.registerPropFunction({'EffectiveElectricalConductivity', fn, {'volumeFraction'} });


            
            fn = @SwellingMaterial.SolidDiffusionModel.updateAverageConcentration;
            % TODO check definition of update function for cAverage in FullSolidDiffusionModel (updateAverageConcentration)
            % If it can be used, then Xavier fix that
            model = model.registerPropFunction({{sd, 'cAverage'}, fn, {{sd, 'c'}}});

            
                      
                   
        end

        
        
        function state = updateRadius(model, state)

            c = state.(sd).cAverage;
            
            radius_0 = model.SolidDiffusion.rp;
            densitySi = model.density;
            MolarMassSi = 28.0855 * 1E-3;
            molarVolumeSi = MolarMassSi/densitySi;
            molarVolumeLi = 8.8 * 1E-6;
            cmaxLi = model.Interface.cmax;
            
            radius = radius_0 * (1 + (3.75*molarVolumeLi*c)/(cmaxLi*molarVolumeSi))^(1/3);
            state.radius = radius;

            if model.use_particle_diffusion
                state.SolidDiffusion.radius = radius;
            end
         end

            
        function state = updatePorosity(model, state, state0, dt)
            porosity = state0.porosity;

            a = state0.volumetricSurfaceArea;       
            R = state0.Interface.R;

            molarVolumeLithiated = model.updateMolarVolumeLithiated(state);
            densitySi = model.Interface.density;
            molarMassSi = 28.0855 * 1E-3;
            molarVolumeSi = molarMassSi/densitySi;

            state.porosity = porosity + dt*a*R*(molarVolumeLithiated - 3.75*molarVolumeSi);
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
            vf = state.Interface.volumeFraction;
            amf = model.activeMaterialFraction;
            radius = state.SolidDiffusion.radius;

            vsa = 3*vf*amf/radius;
            state.Interface.volumetricSurfaceArea = vsa;

            if model.use_particle_diffusion
                state.SolidDiffusion.volumetricSurfaceArea = vsa;
            end
        end

      
         function state = updateMassSource(model, state)
            
            op  = model.operators;
            rp  = state.radius;
            vf  = state.volumeFraction;
            amf = model.activeMaterialFraction;
            
            Rvol = state.SolidDiffusion.Rvol;

            Rvol = op.mapFromBc*Rvol;
            
            state.SolidDiffusion.massSource = - Rvol*((4*pi*rp^3)/(3*amf*vf));
            
         end

         function state = updateEffectiveElectricalConductivity(model, state)
                                 
            sigma = state.Conductivity;
            vf    = state.volumeFraction;
            brugg = model.BruggemanCoefficient;
            
            state.EffectiveElectricalConductivity = sigma*vf.^brugg;
            
        end



        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            time = state0.time + dt;
            state = model.initStateAD(state);

            state = updateControl(model, state, drivingForces);
            

            state                = model.updatePorosity(state, state0, dt);
            state                = model.updateVolumeFraction(state);
            state                = model.updateEffectiveElectricalConductivity(state)         
           

            state                = model.updateStandalonejBcSource(state);
            state                = model.updateCurrent(state);
            state.SolidDiffusion = model.SolidDiffusion.updateMassAccum(state.SolidDiffusion, state0.SolidDiffusion, dt);
            state                = model.dispatchTemperature(state);
            state.SolidDiffusion = model.SolidDiffusion.updateDiffusionCoefficient(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.updateFlux(state.SolidDiffusion);
            state                = model.updateConcentrations(state);
            
            state                = model.updateAverageConcentration(state);
            state                = model.updateRadius(state);
            state                = model.updateVolumetricSurfaceArea(state);
            state.SolidDiffusion = model.SolidDiffusion.updateOperators(state);

            state                = model.updatePhi(state);
            state.Interface      = model.Interface.updateReactionRateCoefficient(state.Interface);
            state.Interface      = model.Interface.updateOCP(state.Interface);
            state.Interface      = model.Interface.updateEta(state.Interface);
            state.Interface      = model.Interface.updateReactionRate(state.Interface);
            state                = model.updateRvol(state);
            state                = model.updateCurrentSource(state);
            state                = model.updateChargeConservation(state);

            state                = model.updateMassSource(state);

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
            
            names = {'chargeCons', ...
                     'massCons', ...
                     'solidDiffusionEq'};
            
            types = {'cell', 'cell', 'cell'};

            primaryVars = model.getPrimaryVariables();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

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
