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
            R = state.R;
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
