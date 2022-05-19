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
            model = model.registerVarNames(varnames);
            
            fn = @SEIActiveMaterial.assembleInterfaceMassCons;
            inputnames = {{itf, 'R'}, {sr, 'R'}, 'R'};
            model = model.registerPropFunction({'interfaceMassCons', fn, inputnames});
            
            fn = @SEIActiveMaterial.updatePotentialDrop;
            inputnames = {{'R'}, {'SolidElectrodeInterface', 'seiwidth'}};
            model = model.registerPropFunction({{itf, 'externalPotentialDrop'}, fn, inputnames});
            model = model.registerPropFunction({{sr, 'externalPotentialDrop'}, fn, inputnames});
            
            fn = @SEIActiveMaterial.dispatchSEIRate;
            inputnames = {{sr, 'R'}};
            model = model.registerPropFunction({{sei, 'R'}, fn, inputnames});
            
            fn = @SEIActiveMaterial.updateSEISurfaceConcentration;
            inputnames = {{sei, 'c'}};
            model = model.registerPropFunction({{sr, 'c'}, fn, inputnames});
        
        
        end
        
        function state = updateSEISurfaceConcentration(model, state)
            
            op = model.(sei).operators;
            
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            c = state.(sei).c;
            
            state.(sr).c = op.mapToIntBc*c;
            
        end
        
        function state = dispatchSEIRate(model, state)

            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            state.(sei).R = state.(sr).R;
            
        end
        
        function state = assembleInterfaceMassCons(model, state)

            itf = 'Interface';
            sr  = 'SideReaction';
            
            Rint = state.(itf).R;
            Rsei = state.(sr).R;          
            R = state.R;
            
            state.interfaceMassCons = R - Rint - Rsei;
            
        end
        
        
        function state = updatePotentialDrop(model, state)
            
            sei = 'SolidElectrodeInterface';
            itf = 'Interface';
            sr  = 'SideReaction';           
            
            kappaSei = model.(sr).conductivity;
            R     = state.R;
            delta = state.(sei).seiwidth;
            
            %% FIXME : check sign
            dphi = - delta.*R/kappasei;
            
            state.(itf).externalPotentialDrop = dphi;
            state.(sr).externalPotentialDrop = dphi;
            
        end
        
        
        function state = dispatchRtotal(model, state)

            itf = 'Interface';
            sr  = 'SideReaction';
            
            state.(itf).Rtotal = state.R;
            state.(sr).Rtotal = state.R;
            
        end

        function getEquations

            itf = 'Interface';
            sd  = 'SolidDiffusion';
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            state = model.dispatchTemperature(state); % model.(itf).T, SolidDiffusion.T, 

            state = model.updateConcentrations(state); % SolidDiffusion.cSurface, model.(itf).cElectrodeSurface, 
            
            state.(itf = state.(itf) = model.(itf).updateOCP(state.(itf)); % model.(itf).OCP, 
            state.(itf) = model.(itf).updateReactionRateCoefficient(state.(itf)); % model.(itf).j0, 

            state = model.updatePotentialDrop(state); % model.(itf).externalPotentialDrop, SideReaction.externalPotentialDrop,
            
            state = model.updatePhi(state); % model.(itf).phiElectrode, 

            state.(itf) = model.(itf).updateEtaWithEx(state.(itf)); % model.(itf).eta, 
            state.(itf) = model.(itf).updateReactionRate(state.(itf)); % model.(itf).R, 
            state = model.dispatchRate(state); % SolidDiffusion.R, 
            state.(sd) = model.(sd).updateMassSource(state.(sd)); % SolidDiffusion.massSource, 
            state.(sd) = model.(sd).assembleSolidDiffusionEquation(state.(sd)); % SolidDiffusion.solidDiffusionEq, 
            
            state = model.updateSEISurfaceConcentration(state);
            state.(sr) = model.(sr).updateReactionRate(state.(sr)); % SideReaction.R, 

            state = model.assembleInterfaceMassCons(state); % interfaceMassCons, 

            state.(sd) = model.(sd).updateDiffusionCoefficient(state.(sd)); % SolidDiffusion.D, 
            state.(sd) = model.(sd).updateFlux(state.(sd)); % SolidDiffusion.flux, 

            % FIXME : add SolidDiffusion.massAccum 

            state.(sd) = model.(sd).updateMassConservation(state.(sd)); % SolidDiffusion.massCons, 

            state = model.updateCurrentSource(state); % eSource, 
            state = model.updateCurrent(state); % j, 
            state = model.updateStandalonejBcSource(state); % jBcSource, 
            state = model.updateChargeConservation(state); % chargeCons, 
            
            state = model.dispatchSEIRate(state); % SolidElectrodemodel.(itf).R, 
            state.(sei) = model.(sei).updateSEIgrowthVelocity(state.(sei)); % SolidElectrodeInterface.v, 
            state.(sei) = model.(sei).updateFlux(state.(sei)); % SolidElectrodeInterface.flux, 

            % FIXME : add sei.massaccum

            state.(sei) = model.(sei).updateMassSource(state.(sei)); % SolidElectrodeInterface.massSource, 
            state.(sei) = model.(sei).updateMassConservation(state.(sei)); % SolidElectrodeInterface.massCons, 

            state.(sei) = model.(sei).assembleWidthEquation(state.(sei), state0.(sei), dt); % SolidElectrodeInterface.widthEq, 
            
            state.(sei) = model.(sei).assembleSolidDiffusionEquation(state.(sei)); % SolidElectrodeInterface.solidDiffusionEq, 
            
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
