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
            model.Control     = CcCvControlModel(paramobj.Control);
            
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
        % Some shorthands used for the sub-models
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


        % Some shorthands used for the sub-models
            an    = 'Anode';
            ct    = 'Cathode';
            itf   = 'Interface';
            elyte = 'Electrolyte';

            phi = state.(elyte).phi;

            state.(an).(itf).phiElectrolyte = phi;
            state.(ct).(itf).phiElectrolyte = phi;
        end

        function state = setupElectrolyteMassCons(model, state)

        % Some shorthands used for the sub-models
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
