classdef SEIActiveMaterial < ActiveMaterial
    
    properties
        
        SolidElectrodeInterface
        SideReaction
        
    end
    
    methods
        
        function model = SEIActiveMaterial(paramobj)

            model = model@ActiveMaterial(paramobj);
            model.SideReaction = SideReaction(paramobj.SideReaction);
            model.SolidElectrodeInterface = SolidElectrodeInterface(paramobj.SolidElectrodeInterface);
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ActiveMaterial(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            varnames = {};
            % Total reaction rate
            varnames{end + 1} = 'R';
            % SEI charge consservation equation
            varnames{end + 1} = 'seiInterfaceChargeCons';
            model = model.registerVarNames(varnames);
            
            fn = @SEIActiveMaterial.assembleSEIchargeCons;
            inputnames = {{itf, 'R'}, {sr, 'R'}, 'R'};
            model = model.registerPropFunction({'seiInterfaceChargeCons', fn, inputnames});
            
            fn = @SEIActiveMaterial.updatePhi;
            model = model.registerPropFunction({{itf, 'phiElectrode'}, fn, {'phi'}});
            model = model.registerPropFunction({{sr, 'phiElectrode'}, fn, {'phi'}});
            
            fn = @SEIActiveMaterial.updatePotentialDrop;
            inputnames = {{'R'}, {'SolidElectrodeInterface', 'delta'}};
            model = model.registerPropFunction({{itf, 'externalPotentialDrop'}, fn, inputnames});
            model = model.registerPropFunction({{sr, 'externalPotentialDrop'}, fn, inputnames});
            
            fn = @SEIActiveMaterial.dispatchSEIRate;
            inputnames = {{sr, 'R'}};
            model = model.registerPropFunction({{sei, 'R'}, fn, inputnames});
            
            fn = @SEIActiveMaterial.updateSEISurfaceConcentration;
            inputnames = {{sei, 'cInterface'}};
            model = model.registerPropFunction({{sr, 'c'}, fn, inputnames});
        
        
        end
        
        function state = updateSEISurfaceConcentration(model, state)

            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            state.(sr).c = state.(sei).cInterface;
            
        end
        
        function state = dispatchSEIRate(model, state)

            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            state.(sei).R = state.(sr).R;
            
        end

        function state = updatePhi(model, state)
            state = updatePhi@ActiveMaterial(model, state);
            state.SideReaction.phiElectrode = state.phi;
        end         
        
        function state = assembleSEIchargeCons(model, state)

            itf = 'Interface';
            sr  = 'SideReaction';
            
            vsa = model.(itf).volumetricSurfaceArea;
            
            Rint = state.(itf).R;
            Rsei = state.(sr).R;
            R = state.R;
            
            state.seiInterfaceChargeCons = R - Rint - Rsei;
            
        end
        
        function state = updatePotentialDrop(model, state)
            
            sei = 'SolidElectrodeInterface';
            itf = 'Interface';
            sr  = 'SideReaction';           
            
            kappaSei = model.(sr).conductivity;
            F = model.(itf).constants.F;
            
            R     = state.R;
            delta = state.(sei).delta;
            
            dphi = F*delta.*R/kappaSei;
            
            state.(itf).externalPotentialDrop = dphi;
            state.(sr).externalPotentialDrop  = dphi;
            
        end
        
        function state = updateControl(model, state, drivingForces, dt)
            
            G = model.G;
            coef = G.cells.volumes;
            coef = coef./(sum(coef));
            
            state.controlCurrentSource = drivingForces.src.*coef;
            
        end
        
        function state = updateStandalonejBcSource(model, state)
            
            state.jBcSource = state.controlCurrentSource;

        end

        function primaryvarnames = getPrimaryVariables(model)
            
            sd  = 'SolidDiffusion';
            sei = 'SolidElectrodeInterface';
            
            primaryvarnames = {{sd, 'c'}           , ...
                               {sd, 'cSurface'}    , ...
                               {'phi'}             , ...
                               {sei, 'c'}          , ...
                               {sei, 'cInterface'} , ...
                               {sei, 'delta'}      , ...
                               {'R'}};
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state, state0)
            
            cleanState = addStaticVariables@ActiveMaterial(model, cleanState, state);
            
            sei = 'SolidElectrodeInterface';
            sr = 'SideReaction';
            
            cleanState.(sei).cExternal     = state.(sei).cExternal;
            cleanState.(sr).phiElectrolyte = state.(sr).phiElectrolyte;
            
        end

        function state = dispatchTemperature(model, state)

            state = dispatchTemperature@ActiveMaterial(model, state);

            sr  = 'SideReaction';
            
            state.(sr).T = state.T;
            
        end
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)

            itf = 'Interface';
            sd  = 'SolidDiffusion';
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            time = state0.time + dt;
            state = model.initStateAD(state);
            
            state = model.dispatchTemperature(state); % model.(itf).T, SolidDiffusion.T

            state = model.updateConcentrations(state); % SolidDiffusion.cSurface, model.(itf).cElectrodeSurface
            
            state.(itf) = model.(itf).updateOCP(state.(itf)); % model.(itf).OCP
            state.(itf) = model.(itf).updateReactionRateCoefficient(state.(itf)); % model.(itf).j0

            state = model.updatePotentialDrop(state); % model.(itf).externalPotentialDrop, SideReaction.externalPotentialDrop
            
            state = model.updatePhi(state); % model.(itf).phiElectrode

            state.(itf) = model.(itf).updateEtaWithEx(state.(itf)); % model.(itf).eta
            state.(itf) = model.(itf).updateReactionRate(state.(itf)); % model.(itf).R
            state = model.dispatchRate(state); % SolidDiffusion.R
            state.(sd) = model.(sd).updateMassSource(state.(sd)); % SolidDiffusion.massSource
            state.(sd) = model.(sd).assembleSolidDiffusionEquation(state.(sd)); % SolidDiffusion.solidDiffusionEq
            
            state = model.updateSEISurfaceConcentration(state);
            state.(sr) = model.(sr).updateReactionRate(state.(sr)); % SideReaction.R

            state = model.assembleSEIchargeCons(state); % SeiIterfaceChargeCons

            state.(sd) = model.(sd).updateDiffusionCoefficient(state.(sd)); % SolidDiffusion.D
            state.(sd) = model.(sd).updateFlux(state.(sd)); % SolidDiffusion.flux

            state.(sd) = model.(sd).updateMassAccum(state.(sd), state0.(sd), dt); % SolidDiffusion.massAccum
            state.(sd) = model.(sd).updateMassConservation(state.(sd)); % SolidDiffusion.massCons

            state = model.updateControl(state, drivingForces, dt);
            state = model.updateCurrentSource(state);       % eSource
            state = model.updateConductivity(state);        % conductivity
            state = model.updateCurrent(state);             % j
            state = model.updateStandalonejBcSource(state); % jBcSource
            state = model.updateChargeConservation(state);  % chargeCons
            
            state = model.dispatchSEIRate(state); % SolidElectrodemodel.(itf).R
            state.(sei) = model.(sei).updateSEIgrowthVelocity(state.(sei)); % SolidElectrodeInterface.v
            state.(sei) = model.(sei).updateFlux(state.(sei)); % SolidElectrodeInterface.flux

            state.(sei) = model.(sei).updateMassAccumTerm(state.(sei), state0.(sei), dt); % SolidElectrodeInterface.accumTerm 
            state.(sei) = model.(sei).updateMassSource(state.(sei)); % SolidElectrodeInterface.massSource
            state.(sei) = model.(sei).updateMassConservation(state.(sei)); % SolidElectrodeInterface.massCons

            state.(sei) = model.(sei).assembleWidthEquation(state.(sei), state0.(sei), dt); % SolidElectrodeInterface.widthEq
            
            state.(sei) = model.(sei).assembleInterfaceBoundaryEquation(state.(sei)); % SolidElectrodeInterface.solidDiffusionEq

            % FIXME : replace the ad-hoc scalings with the correct ones.
            massConsScaling = model.(sd).constants.F; % Faraday constant
            eqs = {};
            eqs{end + 1} = 1e-7*state.chargeCons;
            eqs{end + 1} = 1e9*state.(sd).massCons;
            eqs{end + 1} = 1e5*state.(sd).solidDiffusionEq.*massConsScaling.*model.(itf).G.cells.volumes/dt;
            eqs{end + 1} = 1e10*state.(sei).massCons;
            eqs{end + 1} = state.(sei).widthEq;
            eqs{end + 1} = state.seiInterfaceChargeCons;
            eqs{end + 1} = 1e10*state.(sei).interfaceBoundaryEq;
            
            names = {'chargeCons'              , ...
                     'sd_massCons'             , ...
                     'sd_solidDiffusionEq'     , ...
                     'sei_massCons'            , ...
                     'sei_width'               , ...
                     'sei_interfaceChargeCons', ...
                     'sei_interfaceBoundaryeq'};
            
            types = repmat({'cell'}, 1, numel(names));
            
            primaryVars = model.getPrimaryVariables();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
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
