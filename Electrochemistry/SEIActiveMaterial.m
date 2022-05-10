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
            
            model.dispatchTemperature;                  % model.(itf).T, 
            model.dispatchTemperature;                  % SolidDiffusion.T, 

            model.updateConcentrations;                 % SolidDiffusion.cSurface, 
            model.updateConcentrations;                 % model.(itf).cElectrodeSurface, 
            
            model.(itf).updateOCP;                      % model.(itf).OCP, 
            model.(itf).updateReactionRateCoefficient;  % model.(itf).j0, 

            model.updatePotentialDrop;                  % model.(itf).externalPotentialDrop, 
            model.updatePotentialDrop;                  % SideReaction.externalPotentialDrop, 
            
            model.updatePhi;                            % model.(itf).phiElectrode, 

            model.(itf).updateEtaWithEx;                % model.(itf).eta, 
            model.(itf).updateReactionRate;             % model.(itf).R, 
            model.dispatchRate;                         % SolidDiffusion.R, 
            model.(sd).updateMassSource;                % SolidDiffusion.massSource, 
            model.(sd).assembleSolidDiffusionEquation;  % SolidDiffusion.solidDiffusionEq, 
 
            model.updateSEISurfaceConcentration;
            model.(sr).updateReactionRate;              % SideReaction.R, 

            model.assembleInterfaceMassCons;            % interfaceMassCons, 

            model.(sd).updateDiffusionCoefficient;      % SolidDiffusion.D, 
            model.(sd).updateFlux;                      % SolidDiffusion.flux, 

            % FIXME : add SolidDiffusion.massAccum 

            model.(sd).updateMassConservation;          % SolidDiffusion.massCons, 

            model.updateCurrentSource;                  % eSource, 
            model.updateCurrent;                        % j, 
            model.updateStandalonejBcSource;            % jBcSource, 
            model.updateChargeConservation;             % chargeCons, 
            
            model.dispatchSEIRate;                      % SolidElectrodemodel.(itf).R, 
            model.(sei).updateSEIgrowthVelocity;        % SolidElectrodeInterface.v, 
            model.(sei).updateFlux;                     % SolidElectrodeInterface.flux, 

            % FIXME : add sei.massaccum

            model.(sei).updateMassSource;               % SolidElectrodeInterface.massSource, 
            model.(sei).updateMassConservation;         % SolidElectrodeInterface.massCons, 

            % FIXME : special syntax
            model.(sei).assembleWidthEquation;          % SolidElectrodeInterface.widthEq, 
            
            model.(sei).assembleSolidDiffusionEquation; % SolidElectrodeInterface.solidDiffusionEq, 
            
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
