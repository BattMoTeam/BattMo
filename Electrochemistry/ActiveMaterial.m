classdef ActiveMaterial < ElectronicComponent
    
    properties
        
        %
        % instance of :class:`Interface <Electrochemistry.Electrodes.Interface>`
        %

        Interface

        SolidDiffusion        

        porosity                      % porosity
        volumeFraction                % Volume fraction of the whole material (binder and so on included)
        activeMaterialFraction        % Volume fraction occupied only by the active material
        electricalConductivity        % Electrical conductivity
        InterDiffusionCoefficient     % Inter particle diffusion coefficient parameter (diffusion between the particles)
        density                       %
        thermalConductivity           % Intrinsic Thermal conductivity of the active component
        specificHeatCapacity          % Specific Heat capacity of the active component

        EffectiveDiffusionCoefficient % 

        EffectiveThermalConductivity  % Effective Thermal Conductivity of the active component
        EffectiveVolumetricHeatCapacity % Effective Heat Capacity of the active component

        diffusionModelType
        
        externalCouplingTerm          % only used in no current collector

        BruggemanCoefficient

        use_particle_diffusion
        use_interparticle_diffusion
        standAlone

        isSwellingMaterial  % Boolean equal to 1 if we consider the swelling of the material


                
    end
    
    methods
        
        function model = ActiveMaterial(paramobj)
        %
        % ``paramobj`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %
            model = model@ElectronicComponent(paramobj);

           
            
            fdnames = {'activeMaterialFraction', ...
                       'thermalConductivity'   , ...
                       'electricalConductivity', ...
                       'density'               , ...
                       'specificHeatCapacity'  , ...
                       'externalCouplingTerm'  , ...
                       'diffusionModelType'    , ...
                       'use_thermal'           , ...
                       'BruggemanCoefficient'  , ...
                       'isSwellingMaterial'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % Setup Interface component
            paramobj.Interface.G = model.G;
            
            model.Interface = Interface(paramobj.Interface);

            diffusionModelType = model.diffusionModelType;

            switch model.diffusionModelType
              case 'simple'
                model.SolidDiffusion = SimplifiedSolidDiffusionModel(paramobj.SolidDiffusion);
                model.InterDiffusionCoefficient = paramobj.InterDiffusionCoefficient;
                model.use_particle_diffusion = true;
                model.use_interparticle_diffusion = true;
              case 'full'
                paramobj.SolidDiffusion.np = model.G.cells.num;
                if model.isSwellingMaterial == 0
                    model.SolidDiffusion = FullSolidDiffusionModel(paramobj.SolidDiffusion);
                else
                    model.SolidDiffusion = FullSolidDiffusionSwellingModel(paramobj.SolidDiffusion);
                end
                model.use_particle_diffusion = true;
                model.use_interparticle_diffusion = false;
              case 'interParticleOnly'
                model.InterDiffusionCoefficient = paramobj.InterDiffusionCoefficient;
                model.use_particle_diffusion = false;
                model.use_interparticle_diffusion = true;
              otherwise
                error('Unknown diffusionModelType %s', diffusionModelType);
            end
            
            nc = model.G.cells.num;

            model.volumeFraction = paramobj.volumeFraction*ones(nc, 1);
            model.porosity       = 1 - model.volumeFraction;

            model = model.setupDependentProperties();

            % standAlone=true for standalone simulation of active material. This flag is used only in
            % registerVarAndPropfuncNames (does not impact simulation). Default is false.
            model.standAlone = false;
            
        end
        
        function model = setupDependentProperties(model)           

            model.volumeFraction = 1 - model.porosity;
            vf = model.volumeFraction;
            brugg = model.BruggemanCoefficient;
            
            % setup effective electrical conductivity using Bruggeman approximation
            model.EffectiveElectricalConductivity = model.electricalConductivity.*vf.^brugg;

            
            if model.use_interparticle_diffusion
                
                interDiff = model.InterDiffusionCoefficient;
                amFrac    = model.activeMaterialFraction;
                
                model.EffectiveDiffusionCoefficient = interDiff.*(vf.*amFrac).^brugg;
                
            end

            if model.use_thermal
                % setup effective thermal conductivity
                model.EffectiveThermalConductivity = model.thermalConductivity.*vf.^brugg;
                model.EffectiveVolumetricHeatCapacity = model.specificHeatCapacity.*vf.*model.density;
            end
            
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {'jCoupling', ...
                        'jExternal', ...
                        'SOC'      , ...
                        'Rvol'};
            model = model.registerVarNames(varnames);            

            if model.use_thermal
                varnames = {'jFaceCoupling', ...
                            'jFaceExternal'};
                model = model.registerVarNames(varnames);
            end
            
           
            if strcmp(model.diffusionModelType, 'simple')
                varnames = {'c'        , ...
                            'massCons' , ...
                            'massAccum', ...
                            'massFlux', ...
                            'massSource'};
                model = model.registerVarNames(varnames);
            end
            
            if model.standAlone
                varnames = {'controlCurrentSource'};
                model = model.registerVarNames(varnames);
                varnames = {{itf, 'SOC'}, ...
                            'jCoupling', ...
                            'jExternal'};
                model = model.removeVarNames(varnames);
                varnames = {'T', ...
                            {itf, 'cElectrolyte'},... 
                            {itf, 'phiElectrolyte'}};
                model = model.registerStaticVarNames(varnames);
                
            end
            
            fn = @ActiveMaterial.updateCurrentSource;
            model = model.registerPropFunction({'eSource', fn, {'Rvol'}});
            
            fn = @ActiveMaterial.updatePhi;
            model = model.registerPropFunction({{itf, 'phiElectrode'}, fn, {'phi'}});
            
            fn = @ActiveMaterial.dispatchTemperature;
            model = model.registerPropFunction({{itf, 'T'}, fn, {'T'}});
            if model.use_particle_diffusion
                model = model.registerPropFunction({{sd, 'T'}, fn, {'T'}});
            end

            switch model.diffusionModelType

              case 'simple'

                fn = @ActiveMaterial.updateConcentrations;
                model = model.registerPropFunction({{sd, 'cAverage'}, fn, {'c'}});
                model = model.registerPropFunction({{itf, 'cElectrodeSurface'}, fn, {{sd, 'cSurface'}}});
                
                fn = @ActiveMaterial.updateMassFlux;
                model = model.registerPropFunction({'massFlux', fn, {'c'}});

                fn = @ActiveMaterial.updateMassSource;
                model = model.registerPropFunction({'massSource', fn, {'Rvol'}});

                fn = @ActiveMaterial.updateMassConservation;
                model = model.registerPropFunction({'massCons', fn, {'massAccum', 'massFlux', 'massSource'}});

                fn = @ActiveMaterial.updateRvol;
                model = model.registerPropFunction({'Rvol', fn, {{itf, 'R'}}});
                model = model.registerPropFunction({{sd, 'Rvol'}, fn, {{itf, 'R'}}});

              case 'full'

                fn = @ActiveMaterial.updateConcentrations;
                model = model.registerPropFunction({{itf, 'cElectrodeSurface'}, fn, {{sd, 'cSurface'}}});

                fn = @ActiveMaterial.updateRvol;
                model = model.registerPropFunction({'Rvol', fn, {{itf, 'R'}}});
                model = model.registerPropFunction({{sd, 'Rvol'}, fn, {{itf, 'R'}}});

                % Not used in assembly
                if model.isSwellingMaterial
                    fn = @SwellingMaterial.updateSOC;
                    model = model.registerPropFunction({'SOC', fn, {{sd, 'cAverage'}, 'volumeFraction'}});

                else
                    fn = @ActiveMaterial.updateSOC;
                    model = model.registerPropFunction({'SOC', fn, {{sd, 'cAverage'}}});
                end
                
              case 'interParticleOnly'

                fn = @ActiveMaterial.updateConcentrations;
                model = model.registerPropFunction({{itf, 'cElectrodeSurface'}, fn, {'c'}});
                
                fn = @ActiveMaterial.updateRvol;
                model = model.registerPropFunction({'Rvol', fn, {{itf, 'R'}}});

                fn = @ActiveMaterial.updateMassFlux;
                model = model.registerPropFunction({'massFlux', fn, {'c'}});
                
                fn = @ActiveMaterial.updateMassSource;
                model = model.registerPropFunction({'massSource', fn, {'Rvol'}});
                
                fn = @ActiveMaterial.updateMassConservation;
                model = model.registerPropFunction({'massCons', fn, {'massAccum', 'massFlux', 'massSource'}});

              otherwise
                
                error('diffusionModelType not recognized.');
                
            end

            if model.standAlone
                
                fn = @ActiveMaterial.updateStandalonejBcSource;
                model = model.registerPropFunction({'jBcSource', fn, {'controlCurrentSource'}});

                fn = @ActiveMaterial.updateControl;
                fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
                model = model.registerPropFunction({'controlCurrentSource', fn, {}});
                
            else

                fn = @ActiveMaterial.updatejBcSource;
                model = model.registerPropFunction({'jBcSource', fn, {'jCoupling', 'jExternal'}});

                if model.use_thermal
                    fn = @ActiveMaterial.updatejFaceBc;
                    model = model.registerPropFunction({'jFaceBc', fn, {'jFaceCoupling', 'jFaceExternal'}});
                end
                
                fn = @ActiveMaterial.updatejExternal;
                model = model.registerPropFunction({'jExternal', fn, {}});
                if model.use_thermal
                    model = model.registerPropFunction({'jFaceExternal', fn, {}});
                end
                
                fn = @ActiveMaterial.updatejCoupling;
                model = model.registerPropFunction({'jCoupling', fn, {}});
                if model.use_thermal
                    model = model.registerPropFunction({'jFaceCoupling', fn, {}});
                end
            end
            
            %% Function called to assemble accumulation terms (functions takes in fact as arguments not only state but also state0 and dt)
            if model.use_particle_diffusion & strcmp(model.diffusionModelType, 'simple') | ~model.use_particle_diffusion
                if isSwellingMaterial
                    fn = @SwellingMaterial.assembleAccumTerm;
                    fn = {fn, @(propfunc) PropFunction.accumFuncCallSetupFn(propfunc)};
                    model = model.registerPropFunction({'massAccum', fn, {'c', 'volumeFraction'}});
                else
                    fn = @ActiveMaterial.assembleAccumTerm;
                    fn = {fn, @(propfunc) PropFunction.accumFuncCallSetupFn(propfunc)};
                    model = model.registerPropFunction({'massAccum', fn, {'c'}});
                end
            end

            % we remove this declaration as it is not used in assembly (otherwise it may be computed but not used)
            model = model.removeVarName('SOC');
           
        end
        
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            time = state0.time + dt;
            state = model.initStateAD(state);

            state = updateControl(model, state, drivingForces);
            
            state                = model.updateStandalonejBcSource(state);
            state                = model.updateCurrent(state);
            state.SolidDiffusion = model.SolidDiffusion.updateMassAccum(state.SolidDiffusion, state0.SolidDiffusion, dt);
            state                = model.dispatchTemperature(state);
            state.SolidDiffusion = model.SolidDiffusion.updateDiffusionCoefficient(state.SolidDiffusion);
            state.SolidDiffusion = model.SolidDiffusion.updateFlux(state.SolidDiffusion);
            state                = model.updateConcentrations(state);
            state                = model.updatePhi(state);
            state.Interface      = model.Interface.updateVolumetricSurfaceArea(state.Interface);
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
            rp    = model.(sd).rp;
            vsf   = model.(sd).volumetricSurfaceArea;
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

        function primaryvarnames = getPrimaryVariables(model)
            
            sd = 'SolidDiffusion';
            
            primaryvarnames = {{sd, 'c'}, ...
                               {sd, 'cSurface'}, ...
                               {'phi'}};
            
        end
        
        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end

        function state = updateControl(model, state, drivingForces)
            
            G = model.G;
            coef = G.cells.volumes;
            coef = coef./(sum(coef));
            
            state.controlCurrentSource = drivingForces.src.*coef;
            
        end
        
        
        function cleanState = addStaticVariables(model, cleanState, state, state0)
            
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            
            itf = 'Interface';
            
            cleanState.T = state.T;
            cleanState.(itf).cElectrolyte   = state.(itf).cElectrolyte;
            cleanState.(itf).phiElectrolyte = state.(itf).phiElectrolyte;
            
            sigma = model.electricalConductivity;
            vf    = model.volumeFraction;
            brugg = model.BruggemanCoefficient;
            
            cleanState.conductivity = sigma*vf.^brugg;
            
        end

        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);
            
        end
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            
            [model, state] = prepareTimestep@BaseModel(model, state, state0, dt, drivingForces);
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
        % [state, report] = updateAfterConvergence@ElectronicComponent(model, state0, state, dt, drivingForces);

        % by not calling the parent method, we do not clean the state s
            report = [];
            
        end
         
        function model = validateModel(model, varargin)
            
        end

        %% assembly functions use in this model
         
        function state = updateRvol(model, state)

            vsa = model.Interface.volumetricSurfaceArea;
            
            Rvol = vsa.*state.Interface.R;
            state.Rvol = Rvol;

            if model.use_particle_diffusion
                state.SolidDiffusion.Rvol = Rvol;
            end
            
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
            vf     = model.volumeFraction;
            amFrac = model.activeMaterialFraction;

            c  = state.c;
            c0 = state0.c;

            state.massAccum = vols.*vf.*amFrac.*(c - c0)/dt;
            
        end

        function state = updateMassSource(model, state)
        % used when diffusionModelType == simple
            
            vols = model.G.cells.volumes;
            
            Rvol = state.Rvol;
            
            state.massSource = - Rvol.*vols;
            
        end
        
        
        function state = updateMassConservation(model, state)
        % Used when diffusionModelType == 'simple' or no particle diffusion
            
            flux = state.massFlux;
            source = state.massSource;
            accum = state.massAccum;

            cons = assembleConservationEquation(model, flux, 0, source, accum);
            
            state.massCons = cons;
            
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

            vf       = model.volumeFraction;
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

            vf       = model.volumeFraction;
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