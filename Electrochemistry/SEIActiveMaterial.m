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

            if model.standAlone
                model = model.setupStandAloneModel();
            end
            
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
            % SEI charge consservation equation
            varnames{end + 1} = {itf, 'externalPotentialDrop'};

            model = model.registerVarNames(varnames);

            if model.standAlone
                model = model.registerStaticVarNames({{sei, 'cExternal'}, ...
                                                      {sr, 'phiElectrolyte'}});
            end

            fn = @SEIActiveMaterial.assembleSEIchargeCons;
            inputnames = {{itf, 'R'}, {sr, 'R'}, 'R'};
            model = model.registerPropFunction({'seiInterfaceChargeCons', fn, inputnames});

            fn = @SEIActiveMaterial.dispatchTemperature;
            model = model.registerPropFunction({{sr, 'T'}, fn, {'T'}});
            
            fn = @SEIActiveMaterial.updatePotentialDrop;
            inputnames = {{'R'}, {'SolidElectrodeInterface', 'delta'}};
            model = model.registerPropFunction({{itf, 'externalPotentialDrop'}, fn, inputnames});
            model = model.registerPropFunction({{sr, 'externalPotentialDrop'}, fn, inputnames});

            fn = @SEIActiveMaterial.Interface.updateEtaWithEx;
            varname = VarName({itf}, 'eta');
            inputvarnames = {VarName({itf}, 'phiElectrolyte'), ...
                             VarName({itf}, 'phiElectrode')  , ...
                             VarName({itf}, 'OCP')           , ...
                             VarName({itf}, 'externalPotentialDrop')};
            modelnamespace = {itf};
            propfunction = PropFunction(varname, fn, inputvarnames, modelnamespace);
            model = model.registerPropFunction(propfunction);
            
            fn = @SEIActiveMaterial.dispatchSEIRate;
            inputnames = {{sr, 'R'}};
            model = model.registerPropFunction({{sei, 'R'}, fn, inputnames});
            
            fn = @SEIActiveMaterial.updateSEISurfaceConcentration;
            inputnames = {{sei, 'cInterface'}};
            model = model.registerPropFunction({{sr, 'c'}, fn, inputnames});

            if model.standAlone

                fn = @SEIActiveMaterial.updatePhi;
                model = model.registerPropFunction({{sr, 'phiElectrode'}, fn, {'E'}});

                fn = @SEIActiveMaterial.updateChargeCons;
                inputnames = {'I', ...
                              'R'};
                model = model.registerPropFunction({'chargeCons', fn, inputnames});

            end

        end


        function model = setupStandAloneModel(model)

            if isempty(model.SideReaction)
                return
            end

            model = setupStandAloneModel@ActiveMaterial(model);
            
        end
        
        function state = dispatchTemperature(model, state)

            state = dispatchTemperature@ActiveMaterial(model, state);
            state.SideReaction.T = state.T;
            
        end

        function state = dispatchE(model, state)

            state.phi = state.E;
            
        end

        function state = setupControl(model, state)

            state.R = state.controlCurrentSource;
            
        end

        function state = updateSEISurfaceConcentration(model, state)

            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            state.(sr).c = state.(sei).cInterface;
            
        end

        function state = updateChargeCons(model, state)
            itf = 'Interface';
            F    = model.(itf).constants.F;
            
            I = state.I;
            R = state.R;

            state.chargeCons = I - R*F;
            
        end
        
        function state = dispatchSEIRate(model, state)

            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            state.(sei).R = state.(sr).R;
            
        end

        function state = updatePhi(model, state)

            state = updatePhi@ActiveMaterial(model, state);
            state.SideReaction.phiElectrode = state.E;
            
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
        
        function cleanState = addStaticVariables(model, cleanState, state, state0)
            
            cleanState = addStaticVariables@ActiveMaterial(model, cleanState, state);
            
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            cleanState.(sei).cExternal     = state.(sei).cExternal;
            cleanState.(sr).phiElectrolyte = state.(sr).phiElectrolyte;
            
        end

        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)

            itf = 'Interface';
            sd  = 'SolidDiffusion';
            sei = 'SolidElectrodeInterface';
            sr  = 'SideReaction';
            
            time = state0.time + dt;
            state = model.initStateAD(state);

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end
            
            %% Setup equations and add some scaling
            F   = model.(sd).constants.F;

            scalingcoef = F;
            state.chargeCons = (1/scalingcoef)*state.chargeCons;
            
            rp  = model.(sd).particleRadius;
            vsa = model.(sd).volumetricSurfaceArea;
            
            scalingcoef = vsa*(4*pi*rp^3/3);
            state.(sd).massCons         = (1/scalingcoef)*state.(sd).massCons;
            state.(sd).solidDiffusionEq = (1/scalingcoef)*state.(sd).solidDiffusionEq;

            Mw  = model.(sei).molecularWeight;
            rho = model.(sei).density;

            scalingcoef = 0.5*Mw/rho;
            state.(sei).widthEq = (1/scalingcoef)*state.(sei).widthEq;
            
            deltaref = 1e-2*model.(sd).particleRadius;

            scalingcoef = deltaref;
            state.(sei).interfaceBoundaryEq = (1/scalingcoef)*state.(sei).interfaceBoundaryEq;
            state.(sei).massCons            = (1/scalingcoef)*state.(sei).massCons;
            
            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end            

            names       = model.equationNames;
            types       = model.equationTypes;
            primaryVars = model.primaryVarNames;
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
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
