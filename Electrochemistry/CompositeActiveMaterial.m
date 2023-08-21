classdef CompositeActiveMaterial < ElectronicComponent
    
    properties


        FirstMaterial
        SecondMaterial        

        porosity                      % porosity
        volumeFraction                % Volume fraction of the whole material (binder and so on included)
        activeMaterialFraction        % Volume fraction occupied only by the active material
        electricalConductivity        % Electrical conductivite
        InterDiffusionCoefficient     % Inter particle diffusion coefficient parameter (diffusion between the particles)
        thermalConductivity           % Intrinsic Thermal conductivity of the active component
        heatCapacity                  % Intrinsic Heat capacity of the active component

        EffectiveDiffusionCoefficient % 

        EffectiveThermalConductivity  % Effective Thermal Conductivity of the active component
        EffectiveHeatCapacity         % Effective Heat Capacity of the active component

        
        externalCouplingTerm          % only used in no current collector

        isRoot
        
    end
    
    methods
        
        function model = CompositeActiveMaterial(paramobj)
        %
        % ``paramobj`` is instance of :class:`CompositeActiveMaterialInputParams <Electrochemistry.CompositeActiveMaterialInputParams>`
        %
            model = model@ElectronicComponent(paramobj);
            
            fdnames = {'volumeFraction'        , ...
                       'thermalConductivity'   , ...
                       'electricalConductivity', ...
                       'heatCapacity'          , ...
                       'externalCouplingTerm'  , ...
                       'use_thermal'};
            
            model = dispatchParams(model, paramobj, fdnames);
            

            model.FirstMaterial = ActiveMaterial(paramobj.FirstMaterial);
            model.SecondMaterial = ActiveMaterial(paramobj.SecondMaterial);
            
            nc = model.G.cells.num;

            model.volumeFraction = paramobj.volumeFraction*ones(nc, 1);
            model.porosity       = 1 - model.volumeFraction;
            
            model.isRoot = false;
            
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);
            
            gr = 'FirstMaterial';
            si = 'SecondMaterial';

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            varnames = {'phi'         , ...
                        'jBcSource'   , ...
                        'eSource'     , ...
                        'j'           , ...
                        'chargeCons'  , ...
                        'conductivity', ...
                        'jCoupling'   , ...
                        'jExternal'};

            mvarnames = cellfun(@(varname) {gr, varname}, varnames, 'uniformoutput', false);
            model = model.removeVarNames(mvarnames);
            mvarnames = cellfun(@(varname) {si, varname}, varnames, 'uniformoutput', false);
            model = model.removeVarNames(mvarnames);


            varnames = {'jCoupling', ...
                        'jExternal'};
            
            model = model.registerVarNames(varnames);

            if model.use_thermal
                varnames = {'jFaceCoupling', ...
                            'jFaceExternal'};
                model = model.registerVarNames(varnames);
            end

            if model.isRoot
                varnames = {'controlCurrentSource', ...
                            'cElectrolyte', ...
                            'phiElectrolyte'};
                model = model.registerVarNames(varnames);
                varnames = {'jCoupling', ...
                            'jExternal'};
                model = model.removeVarNames(varnames);
                varnames = {'cElectrolyte',... 
                            'phiElectrolyte', ...
                            'T'};
                model = model.registerStaticVarNames(varnames);
                
            end

            %% register properties

            fn = @CompositeActiveMaterial.updatejBcSource;
            model = model.registerPropFunction({'jBcSource', fn, {'jCoupling', 'jExternal'}});
            
            if model.use_thermal
                fn = @CompositeActiveMaterial.updatejFaceBc;
                model = model.registerPropFunction({'jFaceBc', fn, {'jFaceCoupling', 'jFaceExternal'}});
            end

            fn = @CompositeActiveMaterial.updateConductivity;
            model = model.registerPropFunction({'conductivity', fn, {}});
            
            fn = @CompositeActiveMaterial.updateCurrentSource;
            inputnames = {{si, 'Rvol'}, {gr, 'Rvol'}};
            model = model.registerPropFunction({'eSource', fn, inputnames});
            
            fn = @CompositeActiveMaterial.updatePhi;
            model = model.registerPropFunction({{si, itf, 'phiElectrode'}, fn, {'phi'}});
            model = model.registerPropFunction({{gr, itf, 'phiElectrode'}, fn, {'phi'}});
            
            fn = @CompositeActiveMaterial.dispatchTemperature;
            model = model.registerPropFunction({{si, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{gr, 'T'}, fn, {'T'}});

            if model.isRoot
                
                fn = @CompositeActiveMaterial.updateStandalonejBcSource;
                model = model.registerPropFunction({'jBcSource', fn, {'controlCurrentSource'}});

                fn = @CompositeActiveMaterial.updateStandAloneElectrolyte;
                inputnames = {'cElectrolyte', 'phiElectrolyte'};
                model = model.registerPropFunction({{si, itf, 'cElectrolyte'}, fn, inputnames});
                model = model.registerPropFunction({{si, itf, 'phiElectrolyte'}, fn, inputnames});
                model = model.registerPropFunction({{gr, itf, 'cElectrolyte'}, fn, inputnames});
                model = model.registerPropFunction({{gr, itf, 'phiElectrolyte'}, fn, inputnames});
                
            else
                
                fn = @CompositeActiveMaterial.updatejExternal;
                model = model.registerPropFunction({'jExternal', fn, {}});
                if model.use_thermal
                    model = model.registerPropFunction({'jFaceExternal', fn, {}});
                end

                fn = @CompositeActiveMaterial.updatejCoupling;
                model = model.registerPropFunction({'jCoupling', fn, {}});
                if model.use_thermal
                    model = model.registerPropFunction({'jFaceCoupling', fn, {}});
                end

            end
            
           
        end

        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)

            gr = 'FirstMaterial';
            si = 'SecondMaterial';
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            time = state0.time + dt;
            state = model.initStateAD(state);

            state = updateControl(model, state, drivingForces, dt);
            
            state = model.updateConductivity(state);
            state = model.updateStandalonejBcSource(state);
            state = model.updateCurrent(state);
            state = model.dispatchTemperature(state);
            state.(si).(sd) = model.(si).(sd).updateMassAccum(state.(si).(sd), state0.(si).(sd), dt);
            state.(si).(sd) = model.(si).(sd).updateAverageConcentration(state.(si).(sd));
            state.(si) = model.(si).updateSOC(state.(si));
            state.(si) = model.(si).dispatchTemperature(state.(si));
            state.(si).(sd) = model.(si).(sd).updateDiffusionCoefficient(state.(si).(sd));
            state.(si).(sd) = model.(si).(sd).updateFlux(state.(si).(sd));
            state = model.updateStandAloneElectrolyte(state);
            state.(si) = model.(si).updateConcentrations(state.(si));
            state = model.updatePhi(state);
            state.(si).(itf) = model.(si).(itf).updateReactionRateCoefficient(state.(si).(itf));
            state.(si).(itf) = model.(si).(itf).updateOCP(state.(si).(itf));
            state.(si).(itf) = model.(si).(itf).updateEta(state.(si).(itf));
            state.(si).(itf) = model.(si).(itf).updateReactionRate(state.(si).(itf));
            state.(si) = model.(si).updateRvol(state.(si));
            state.(si).(sd) = model.(si).(sd).updateMassSource(state.(si).(sd));
            state.(si).(sd) = model.(si).(sd).assembleSolidDiffusionEquation(state.(si).(sd));
            state.(si).(sd) = model.(si).(sd).updateMassConservation(state.(si).(sd));
            state.(gr).(sd) = model.(gr).(sd).updateMassAccum(state.(gr).(sd), state0.(gr).(sd), dt);
            state.(gr).(sd) = model.(gr).(sd).updateAverageConcentration(state.(gr).(sd));
            state.(gr) = model.(gr).updateSOC(state.(gr));
            state.(gr) = model.(gr).dispatchTemperature(state.(gr));
            state.(gr).(sd) = model.(gr).(sd).updateDiffusionCoefficient(state.(gr).(sd));
            state.(gr).(sd) = model.(gr).(sd).updateFlux(state.(gr).(sd));
            state.(gr) = model.(gr).updateConcentrations(state.(gr));
            state.(gr).(itf) = model.(gr).(itf).updateReactionRateCoefficient(state.(gr).(itf));
            state.(gr).(itf) = model.(gr).(itf).updateOCP(state.(gr).(itf));
            state.(gr).(itf) = model.(gr).(itf).updateEta(state.(gr).(itf));
            state.(gr).(itf) = model.(gr).(itf).updateReactionRate(state.(gr).(itf));
            state.(gr) = model.(gr).updateRvol(state.(gr));
            state = model.updateCurrentSource(state);
            state = model.updateChargeConservation(state);
            state.(gr).(sd) = model.(gr).(sd).updateMassSource(state.(gr).(sd));
            state.(gr).(sd) = model.(gr).(sd).assembleSolidDiffusionEquation(state.(gr).(sd));
            state.(gr).(sd) = model.(gr).(sd).updateMassConservation(state.(gr).(sd));
            
            %% Setup equations and add some scaling

            
            eqs = {};
            eqs{end + 1} = state.chargeCons;

            % Equations for graphite
            n     = model.(gr).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F     = model.(gr).(sd).constants.F;
            vol   = model.operators.pv;
            rp    = model.(gr).(sd).rp;
            vsf   = model.(gr).(sd).volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
            eqs{end + 1} = scalingcoef*state.(gr).(sd).massCons;
            eqs{end + 1} = scalingcoef*state.(gr).(sd).solidDiffusionEq;

            % Equations for SecondMaterial
            n     = model.(si).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F     = model.(si).(sd).constants.F;
            vol   = model.operators.pv;
            rp    = model.(si).(sd).rp;
            vsf   = model.(si).(sd).volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
            eqs{end + 1} = scalingcoef*state.(si).(sd).massCons;
            eqs{end + 1} = scalingcoef*state.(si).(sd).solidDiffusionEq;

            
            names = {'chargeCons'         , ...
                     'gr_massCons'        , ...
                     'gr_solidDiffusionEq', ...
                     'si_massCons'        , ...
                     'si_solidDiffusionEq'};
            
            types = {'cell', 'cell', 'cell', 'cell', 'cell'};

            primaryVars = model.getPrimaryVariableNames();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function primaryvarnames = getPrimaryVariableNames(model)
            
            gr = 'FirstMaterial';
            si = 'SecondMaterial';

            sd = 'SolidDiffusion';

            primaryvarnames = {{gr, sd, 'c'}       , ...
                               {gr, sd, 'cSurface'}, ...
                               {si, sd, 'c'}       , ...
                               {si, sd, 'cSurface'}, ...
                               {'phi'}};
            
        end
        
        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end

        function state = updateControl(model, state, drivingForces, dt)
            
            op = model.op;
            coef = op.getCellVolumes();
            coef = coef./(sum(coef));
            
            state.controlCurrentSource = drivingForces.src.*coef;
            
        end

        function cleanState = addStaticVariables(model, cleanState, state, state0)
            
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            
            cleanState.T = state.T;
            cleanState.cElectrolyte   = state.cElectrolyte;
            cleanState.phiElectrolyte = state.phiElectrolyte;
            
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

        function state = updateConductivity(model, state)

            G  = model.G;
            vf = model.volumeFraction;
            
            gr  = 'FirstMaterial';
            si  = 'SecondMaterial';

            nc = G.cells.num;

            mats = {gr, si};

            sigma = zeros(nc, 1);
            for imat = 1 : numel(mats)
                mat = mats{imat};
                amvf     = model.(mat).activeMaterialFraction;
                bg       = model.(mat).BruggemanCoefficient;
                sigmamat = model.(mat).electricalConductivity;
                sigma = ((vf.*amvf).^bg).*sigmamat;
            end

            state.conductivity = sigma;
            
        end
        
        
        function state = updateStandalonejBcSource(model, state)
            
            state.jBcSource = state.controlCurrentSource;

        end

        function state = updateCurrentSource(model, state)

            gr  = 'FirstMaterial';
            si  = 'SecondMaterial';
            itf = 'Interface';
            
            %% TODO : check whether constants n should always be the same for graphite and silicon (and impose them from parent)
            vols = model.operators.getCellVolumes();
            F    = model.(gr).(itf).constants.F;
            nGr  = model.(gr).(itf).n;
            nSi  = model.(si).(itf).n;
            
            state.eSource = - (nGr*state.(gr).Rvol + nSi*state.(si).Rvol).*vols.*F;
            
        end
        
        
        function state = updatePhi(model, state)
            
            gr = 'FirstMaterial';
            si = 'SecondMaterial';
            
            itf = 'Interface';

            state.(gr).(itf).phiElectrode = state.phi;
            state.(si).(itf).phiElectrode = state.phi;
            
        end         
        
        function state = dispatchTemperature(model, state)

            gr = 'FirstMaterial';
            si = 'SecondMaterial';
            
            state.(gr).T = state.T;
            state.(si).T = state.T;
            
        end


        function state = updateStandAloneElectrolyte(model, state)

            gr = 'FirstMaterial';
            si = 'SecondMaterial';
            
            itf = 'Interface';

            state.(gr).(itf).cElectrolyte   = state.cElectrolyte;
            state.(gr).(itf).phiElectrolyte = state.phiElectrolyte;
            state.(si).(itf).cElectrolyte   = state.cElectrolyte;
            state.(si).(itf).phiElectrolyte = state.phiElectrolyte;

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
        

        function stop = stopfunction(model, state, state0_inner)
            
            gr = 'FirstMaterial';
            si = 'SecondMaterial';

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            mats = {gr, si};
            
            stop = false;
            
            for imat = 1 : numel(mats)
                mat = mats{imat};
                cmin = (model.(mat).(itf).theta0)*(model.(mat).(itf).cmax);
                vols = model.(mat).(sd).operators.vols;

                % In following function, we assume that we have only one particle
                c = state.(mat).(sd).c;
                cAverage = (sum(vols.*c)/sum(vols));

                if cAverage <= cmin
                    stop = true;
                    return
                end

            end

        end


    end
    
end


%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
