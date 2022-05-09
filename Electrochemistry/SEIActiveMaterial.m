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
            
            SolidElectrodeInterface.updateFlux,       SolidElectrodeInterface.flux, 
            SolidElectrodeInterface.updateSEIgrowthVelocity,       SolidElectrodeInterface.v, 
            SolidElectrodeInterface.updateMassConservation,       SolidElectrodeInterface.massCons, 
            SolidElectrodeInterface.updateMassSource,       SolidElectrodeInterface.massSource, 
            SolidElectrodeInterface.assembleSolidDiffusionEquation,       SolidElectrodeInterface.solidDiffusionEq, 
            SolidElectrodeInterface.assembleWidthEquation,       SolidElectrodeInterface.widthEq, 
            SideReaction.updateReactionRate,       SideReaction.R, 
            Interface.updateReactionRateCoefficient,       Interface.j0, 
            Interface.updateOCP,       Interface.OCP, 
            Interface.updateEtaWithEx,       Interface.eta, 
            Interface.updateReactionRate,       Interface.R, 
            SolidDiffusionModel.updateDiffusionCoefficient,       SolidDiffusion.D, 
            SolidDiffusionModel.updateFlux,       SolidDiffusion.flux, 
            SolidDiffusionModel.updateMassConservation,       SolidDiffusion.massCons, 
            SolidDiffusionModel.updateMassSource,       SolidDiffusion.massSource, 
            SolidDiffusionModel.assembleSolidDiffusionEquation,       SolidDiffusion.solidDiffusionEq, 
            ElectronicComponent.updateCurrent,       j, 
            ElectronicComponent.updateChargeConservation,       chargeCons, 
            ActiveMaterial.updateStandalonejBcSource,       jBcSource, 
            ActiveMaterial.updateCurrentSource,       eSource, 
            ActiveMaterial.updatePhi,       Interface.phiElectrode, 
            ActiveMaterial.dispatchTemperature,       Interface.T, 
            ActiveMaterial.dispatchTemperature,       SolidDiffusion.T, 
            ActiveMaterial.updateConcentrations,       SolidDiffusion.cSurface, 
            ActiveMaterial.updateConcentrations,       Interface.cElectrodeSurface, 
            ActiveMaterial.dispatchRate,       SolidDiffusion.R, 
            SEIActiveMaterial.assembleInterfaceMassCons,       interfaceMassCons, 
            SEIActiveMaterial.updatePotentialDrop,       Interface.externalPotentialDrop, 
            SEIActiveMaterial.updatePotentialDrop,       SideReaction.externalPotentialDrop, 
            SEIActiveMaterial.dispatchSEIRate,       SolidElectrodeInterface.R, 
            
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
