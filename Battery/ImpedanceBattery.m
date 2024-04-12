classdef ImpedanceBattery < GenericBattery

    properties

    end

    methods

        function model = ImpedanceBattery(inputparams)

            model = model@GenericBattery(inputparams);

        end


        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@GenericBattery(model);

            % defines shorthands for the submodels
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            am1     = 'ActiveMaterial1';
            am2     = 'ActiveMaterial2';
            cc      = 'CurrentCollector';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            thermal = 'ThermalModel';
            ctrl    = 'Control';
            sr      = 'SideReaction';
            sei     = 'SolidElectrodeInterface';

            model = model.removeVarName({ctrl, 'ctrlType'});
            
            fn = @ImpedendanceBattery.updateNegativeElectrodeMassAccum;
            inputvarnames = {{ne, co, am, sd, 'c'}, {ctrl, 'omega'}};
            outputvarname = {ne, co, am, sd, 'massAccum'};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ImpedendanceBattery.updatePositiveElectrodeMassAccum;
            inputvarnames = {{pe, co, am, sd, 'c'}, {ctrl, 'omega'}};
            outputvarname = {pe, co, am, sd, 'massAccum'};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ImpedendanceBattery.updateElectrolyteMassAccum;
            inputvarnames = {{elyte, 'c'}, {ctrl, 'omega'}};
            outputvarname = {elyte, 'massAccum'};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
        end

        function Z = computeImpedance(model, state, omega)

            ctrl = 'Control';

            state.(ctrl).omega = omega;
            state.(ctrl).I = 0;
            
            state = model.initStateAD(state);

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            state = model.applyScaling(state);
            
            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end
            eqs = vertcat(eqs{:});
            
            indI = model.getIndexPrimaryVariable({ctrl, 'I'});

            inds = true(numel(model.getPrimaryVariableNames), 1);
            inds(indI) = false;
            inds = find(inds);

            b = eqs.jac{indI};

            eqs.jac = eqs.jac(inds); 
            eqs = combineEquations(eqs);

            inds = model.getRangePrimaryVariable(state.(ctrl).E, {ctrl, 'E'});
            
            A = eqs.jac{1};
            
            x = A\b;
            
            Z = x(inds(1) : inds(2));
            
        end


        function state = updateNegativeElectrodeMassAccum(model, state)
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';

            op = model.(ne).(co).(am).(sd).operators;

            c     = state.(ne).(co).(am).(sd).c;
            omega = state.Control.omega;
            
            state.(ne).(co).(am).(sd).massAccum = i*omega.*op.vols.*c;

            
        end
        
        function state = updatePositiveElectrodeMassAccum(model, state)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';

            op = model.(pe).(co).(am).(sd).operators;

            c     = state.(pe).(co).(am).(sd).c;
            omega = state.Control.omega;
            
            state.(pe).(co).(am).(sd).massAccum = i*omega.*op.vols.*c;

        end
        
        function state = updateElectrolyteMassAccum(model, state)

            elyte   = 'Electrolyte';
            
            c     = state.(elyte).c;
            omega = state.Control.omega;

            effectiveVolumes = model.(elyte).volumeFraction.*model.(elyte).G.getVolumes();

            state.(elyte).massAccum  = i*omega.*effectiveVolumes.*c;
            
        end                
    end

    methods (Static)


        function isequal = compareVarName(varname1, varname2)
        % This generic utility function could have been moved another place than in the impedance model.
            
            if isempty(varname1) & isempty(varname2)
                % Needed in recursion
                isequal = true;
                return
            end
            
            if numel(varname1) == numel(varname2)

                svarname1 = varname1{1};
                svarname2 = varname2{1};

                if ischar(svarname1) & ischar(svarname2)
                    if strcmp(svarname1, svarname2)
                        isequal = ImpedanceBattery.compareVarName(varname1(2 : end), varname2(2 : end));
                    else
                        isequal = false;
                        return
                    end
                elseif isnumeric(svarname1) & isnumeric(svarname2)
                    if svarname1 == svarname2
                        isequal = ImpedanceBattery.compareVarName(varname1(2 : end), varname2(2 : end));
                    else
                        isequal = false;
                        return
                    end
                else
                    isequal = false;
                    return
                end
            else
                isequal = false;
                return
            end
        end
    end

end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
