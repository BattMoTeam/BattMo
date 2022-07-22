classdef SingleParticleSEI < BaseModel
    
    properties
        
        Anode
        Cathode
        Electrolyte
        Control
        
    end
    
    methods
        
        function model = SingleParticleSEI(paramobj)

            model = model@BaseModel();

            model.Anode       = SEIActiveMaterial(paramobj.Anode);
            model.Cathode     = ActiveMaterial(paramobj.Cathode);
            model.Electrolyte = SingleCellElectrolyte(paramobj.Electrolyte);
            % model.Control    = CcCvControlModel(paramobj.Control);
            model.Control     = ControlModel(paramobj.Control);
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % temperature
            varnames{end + 1} = 'T';
            model = model.registerVarNames(varnames);
            
            % Some shorthands used for the sub-models
            an    = 'Anode';
            ct    = 'Cathode';
            sd    = 'SolidDiffusion';
            itf   = 'Interface';
            sei   = 'SolidElectrodeInterface';
            sr    = 'SideReaction';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
            
            fn = @SingleParticleSEI.dispatchEcathode;
            model = model.registerPropFunction({{ct, 'phi'}, fn, {{ctrl, 'E'}}});

            fn = @SingleParticleSEI.dispatchT;
            model = model.registerPropFunction({{an, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{ct, 'T'}, fn, {'T'}});
            
            fn = @SingleParticleSEI.setupEIequation;
            model = model.registerPropFunction({{ctrl, 'I'}, fn, {{an, itf, 'R'}}});

            fn = @SingleParticleSEI.setupElectrolyteCoupling;
            model = model.registerPropFunction({{an, itf, 'phiElectrolyte'}, fn, {{elyte, 'phi'}}});
            model = model.registerPropFunction({{an, sr, 'phiElectrolyte'}, fn, {{elyte, 'phi'}}});
            model = model.registerPropFunction({{ct, itf, 'phiElectrolyte'}, fn, {{elyte, 'phi'}}});
            
            fn = @SingleParticleSEI.setupElectrolyteMassCons;
            model = model.registerPropFunction({{elyte, 'massCons'}, fn, {{an, 'R'}, {ct, itf, 'R'}}});

            eldes = {an, ct};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                switch elde
                  case an
                    fn = @SingleParticleSEI.updateAnodeReactionRateCoefficient;
                  case ct
                    fn = @SingleParticleSEI.updateCathodeReactionRateCoefficient;
                  otherwise
                    error('elde not recognized');
                end
                inputnames = {VarName({elde, itf}, 'cElectrodeSurface'), ...
                              VarName({elde, itf}, 'T')};
                modelnamespace = {elde, itf};
                varname = VarName({elde, itf}, 'j0');
                propfunction = PropFunction(varname, fn, inputnames, modelnamespace);
                model = model.registerPropFunction(propfunction);
            end


            varnames = {'chargeCons', 'j', 'jBcSource', 'jCoupling', 'jBcSource', 'conductivity', 'eSource'};
            for ivar = 1 : numel(varnames)
                model = model.removeVarName({ct, varnames{ivar}});
                model = model.removeVarName({an, varnames{ivar}});
            end

            varnames = {'dUdT', 'cElectrolyte', 'SOC'};
            for ivar = 1 : numel(varnames)
                model = model.removeVarName({ct, itf, varnames{ivar}});
                model = model.removeVarName({an, itf, varnames{ivar}});
            end
            
        end

        function state = dispatchT(model, state)

            ct = 'Cathode';
            an = 'Anode';

            state.(ct).T = state.T;
            state.(an).T = state.T;
            
        end
        
        function state = dispatchEcathode(model, state)

            ct   = 'Cathode';
            ctrl = 'Control';
            
            E = state.(ctrl).E

            state.(ct).phi = E;
            
        end

        function state = setupEIequation(model, state)

            an   = 'Anode';
            itf  = 'Interface';
            ctrl = 'Control';

            F = model.(an).constants.F;
            anArea = model.anodeArea; 
            
            R = state.(an).(itf).R;
            I = state.(ctrl).I;

            staet.(ctrl).EIequation = I - anArea*F*R;
            
        end

        function state = setupElectrolyteCoupling(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            itf   = 'Interface';
            elyte = 'Electrolyte';

            phi = state.(elyte).phi;

            state.(an).(itf).phiElectrolyte = phi;
            state.(ct).(itf).phiElectrolyte = phi;
        end

        function state = setupElectrolyteMassCons(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            itf   = 'Interface';
            elyte = 'Electrolyte';

            anArea = model.anodeArea;
            ctArea = model.cathodeArea;

            anR = state.(an).R;
            ctR = state.(ct).(itf).R;
            
            state.(elyte).massCons = anR*anArea + ctR*ctArea;
            
        end
        

        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)

            an    = 'Anode';
            ct    = 'Cathode';
            sd    = 'SolidDiffusion';
            itf   = 'Interface';
            sei   = 'SolidElectrodeInterface';
            sr    = 'SideReaction';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            time = state0.time + dt;
            state = model.initStateAD(state);
            
            state            = model.dispatchEcathode(state);
            state            = model.dispatchT(state);
            state.(ct)       = model.(ct).dispatchTemperature(state.(ct));
            state.(ct).(sd)  = model.(ct).(sd).updateDiffusionCoefficient(state.(ct).(sd));
            state.(ct).(sd)  = model.(ct).(sd).updateFlux(state.(ct).(sd));
            state.(ct)       = model.(ct).updateConcentrations(state.(ct));
            state            = model.setupElectrolyteCoupling(state);
            state.(ct)       = model.(ct).updateConcentrations(state.(ct));
            state.(ct)       = model.(ct).updatePhi(state.(ct));
            state.(ct)       = model.(ct).dispatchTemperature(state.(ct));
            state.(ct).(itf) = model.(ct).(itf).updateCathodeReactionRateCoefficient(state.(ct).(itf));
            state.(ct).(itf) = model.(ct).(itf).updateOCP(state.(ct).(itf));
            state.(ct).(itf) = model.(ct).(itf).updateEtaWithEx(state.(ct).(itf));
            state.(ct).(itf) = model.(ct).(itf).updateReactionRate(state.(ct).(itf));
            state.(ct)       = model.(ct).dispatchRate(state.(ct));
            state.(ct).(sd)  = model.(ct).(sd).updateMassSource(state.(ct).(sd));
            state.(ct).(sd)  = model.(ct).(sd).assembleSolidDiffusionEquation(state.(ct).(sd));
            state.(ct).(sd)  = model.(ct).(sd).updateMassConservation(state.(ct).(sd));
            state            = model.setupElectrolyteMassCons(state);
            state            = model.dispatchT(state);
            state.(an)       = model.(an).dispatchTemperature(state.(an));
            state.(an).(sd)  = model.(an).(sd).updateDiffusionCoefficient(state.(an).(sd));
            state.(an).(sd)  = model.(an).(sd).updateFlux(state.(an).(sd));
            state.(an)       = model.(an).updateConcentrations(state.(an));
            state            = model.setupElectrolyteCoupling(state);
            state.(an)       = model.(an).updateConcentrations(state.(an));
            state.(an)       = model.(an).updatePhi(state.(an));
            state.(an)       = model.(an).dispatchTemperature(state.(an));
            state.(an).(itf) = model.(an).(itf).updateAnodeReactionRateCoefficient(state.(an).(itf));
            state.(an).(itf) = model.(an).(itf).updateOCP(state.(an).(itf));
            state            = model.setupElectrolyteCoupling(state);
            state.(an)       = model.(an).updatePhi(state.(an));
            state.(an)       = model.(an).updatePotentialDrop(state.(an));
            state.(an).(itf) = model.(an).(itf).updateEtaWithEx(state.(an).(itf));
            state.(an).(itf) = model.(an).(itf).updateReactionRate(state.(an).(itf));
            state            = model.setupEIequation(state);
            state.(ctrl)     = model.(ctrl).updateControlEquation(state.(ctrl));
            state.(an)       = model.(an).dispatchRate(state.(an));
            state.(an).(sd)  = model.(an).(sd).updateMassSource(state.(an).(sd));
            state.(an).(sd)  = model.(an).(sd).assembleSolidDiffusionEquation(state.(an).(sd));
            state.(an).(sd)  = model.(an).(sd).updateMassConservation(state.(an).(sd));
            state.(an)       = model.(an).updatePotentialDrop(state.(an));
            state.(an)       = model.(an).updateSEISurfaceConcentration(state.(an));
            state.(an).(sr)  = model.(an).(sr).updateReactionRate(state.(an).(sr));
            state.(an)       = model.(an).assembleSEIchargeCons(state.(an));
            state.(an)       = model.(an).dispatchSEIRate(state.(an));
            state.(an).(sei) = model.(an).(sei).updateSEIgrowthVelocity(state.(an).(sei));
            state.(an).(sei) = model.(an).(sei).assembleWidthEquation(state.(an).(sei));
            state.(an).(sei) = model.(an).(sei).updateMassSource(state.(an).(sei));
            state.(an).(sei) = model.(an).(sei).assembleInterfaceBoundaryEquation(state.(an).(sei));
            state.(an).(sei) = model.(an).(sei).updateFlux(state.(an).(sei));
            state.(an).(sei) = model.(an).(sei).updateMassAccumTerm(state.(an).(sei));
            state.(an).(sei) = model.(an).(sei).updateMassConservation(state.(an).(sei));

            eqs = {};
            eqs{end + 1} = state.(an).(sd).massCons;
            eqs{end + 1} = state.(an).(sd).solidDiffusionEq;
            eqs{end + 1} = state.(an).(sei).massCons;
            eqs{end + 1} = state.(an).(sei).widthEq;
            eqs{end + 1} = state.(an).seiInterfaceChargeCons;
            eqs{end + 1} = state.(an).(sei).interfaceBoundaryEq;
            eqs{end + 1} = state.(ct).(sd).massCons;
            eqs{end + 1} = state.(ct).(sd).solidDiffusionEq;
            eqs{end + 1} = state.(elyte).massCons;
            eqs{end + 1} = state.(ctrl).controlEquation;
            
            names = {'an_sd_massCons'             , ...
                     'an_sd_solidDiffusionEq'     , ...
                     'an_sei_massCons'            , ...
                     'an_sei_widthEq'             , ...
                     'an_seiInterfaceChargeCons' , ...
                     'an_sei_interfaceBoundaryEq' , ...
                     'ct_sd_massCons'             , ...
                     'ct_sd_solidDiffusionEq'     , ...
                     'elyte_massCons'             , ...
                     'ctrl_controlEq'};
            
            types = repmat({'cell'}, 1, numel(names));
            
            primaryVars = model.getPrimaryVariables();
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end

        function primaryvarnames = getPrimaryVariables(model)
            
            an    = 'Anode';
            ct    = 'Cathode';
            sd    = 'SolidDiffusion';
            itf   = 'Interface';
            sei   = 'SolidElectrodeInterface';
            sr    = 'SideReaction';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
            
            primaryvarnames = {{an, sd, 'c'}           , ...
                               {an, sd, 'cSurface'}    , ...
                               {an, sei, 'c'}          , ...
                               {an, sei, 'cInterface'} , ...
                               {an, sei, 'delta'}      , ...
                               {ct, sd, 'c'}           , ...
                               {ct, sd, 'cSurface'}    , ...
                               {elyte, 'phi'}          , ...
                               {ctrl, 'I'}             , ...
                               {ctrl, 'E'}};
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state, state0)

            cleanState = addStaticVariables@BaseModel(model, cleanState, state);

            an  = 'Anode';
            ct  = 'Cathode';
            sei = 'SolidElectrodeInterface';
            itf = 'Interface';
            sr  = 'SideReaction';

            cleanState.T = state.T;
            % boundary condition
            cleanState.(an).phi = 0;
            % no potential drop at the cathode
            cleanState.(ct).(itf).externalPotentialDrop = 0;
            % external EC concentration is set to constant
            cleanState.(an).(sei).cExternal = state.(an).(sei).cExternal;
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
