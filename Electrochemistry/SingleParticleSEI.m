classdef SingleParticleSEI < BaseModel
    
    properties
        
        Anode
        Cathode
        Electrolyte
        Control
   
        anodeArea   % anode area / [m^2]
        cathodeArea % cathode area / [m^2]

    end
    
    methods
        
        function model = SingleParticleSEI(inputparams)

            model = model@BaseModel();

            model.Anode       = SEIActiveMaterial(inputparams.Anode);
            model.Cathode     = ActiveMaterial(inputparams.Cathode);
            model.Electrolyte = SingleCellElectrolyte(inputparams.Electrolyte);
            model.Control     = model.setupControl(inputparams.Control);
            
            fdnames = {'anodeArea', ...
                       'cathodeArea'};
            
            model = dispatchParams(model, inputparams, fdnames);
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            % Some shorthands used for the sub-models
            an    = 'Anode';
            ct    = 'Cathode';
            sd    = 'SolidDiffusion';
            itf   = 'Interface';
            sei   = 'SolidElectrodeInterface';
            sr    = 'SideReaction';
            elyte = 'Electrolyte';
            ctrl  = 'Control';

            varnames = {};
            % temperature
            varnames{end + 1} = 'T';
            % time
            varnames{end + 1} = 'time';
            % potential at the cathode
            varnames{end + 1} = {ct, 'phi'};
            model = model.registerVarNames(varnames);
            
            fn = @SingleParticleSEI.dispatchEcathode;
            model = model.registerPropFunction({{ct, 'phi'}, fn, {{ctrl, 'E'}}});

            fn = @SingleParticleSEI.dispatchT;
            model = model.registerPropFunction({{an, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{ct, 'T'}, fn, {'T'}});
            
            fn = @SingleParticleSEI.dispatchTime;
            model = model.registerPropFunction({{ctrl, 'time'}, fn, {'time'}});
            
            fn = @SingleParticleSEI.setupEIequation;
            model = model.registerPropFunction({{ctrl, 'EIequation'}, fn, {{an, itf, 'intercalationFlux'}, {ctrl, 'I'}}});

            fn = @SingleParticleSEI.setupElectrolyteCoupling;
            model = model.registerPropFunction({{an, itf, 'phiElectrolyte'}, fn, {{elyte, 'phi'}}});
            model = model.registerPropFunction({{an, sr, 'phiElectrolyte'}, fn, {{elyte, 'phi'}}});
            model = model.registerPropFunction({{ct, itf, 'phiElectrolyte'}, fn, {{elyte, 'phi'}}});
            
            fn = @SingleParticleSEI.setupElectrolyteMassCons;
            model = model.registerPropFunction({{elyte, 'massCons'}, fn, {{an, 'intercalationFlux'}, {ct, itf, 'intercalationFlux'}}});

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
                modelnamespace = {};
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
        
        function control = setupControl(model, inputparams)

            switch inputparams.controlPolicy
              case "CCDischarge"
                control = CCDischargeControlModel(inputparams); 
              case "CCCV"
                control = seiCcCvControlModel(inputparams);
              case "CV"
                control = CvControlModel(inputparams);
              otherwise
                error('Error controlPolicy not recognized');
            end
            
        end


        function state = dispatchTime(model, state)

            state.Control.time = state.time;
            
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
            
            E = state.(ctrl).E;

            state.(ct).phi = E;
            
        end

        function state = setupEIequation(model, state)

            an   = 'Anode';
            itf  = 'Interface';
            ctrl = 'Control';

            F = model.(an).constants.F;
            anArea = model.anodeArea; 
           
            R = state.(an).R;
            I = state.(ctrl).I;
            
            state.(ctrl).EIequation = I -  anArea*F*R;
            
        end

        function state = setupElectrolyteCoupling(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            itf   = 'Interface';
            sr    = 'SideReaction';
            elyte = 'Electrolyte';

            phi = state.(elyte).phi;

            state.(an).(itf).phiElectrolyte = phi;
            state.(an).(sr).phiElectrolyte  = phi;
            state.(ct).(itf).phiElectrolyte = phi;
        end

        function state = updateAnodeReactionRateCoefficient(model, state)

            an  = 'Anode';
            itf = 'Interface';

            state.(an).(itf) = model.updateReactionRateCoefficient(model.(an).(itf), state.(an).(itf));
            
        end

        function state = updateCathodeReactionRateCoefficient(model, state)

            ct  = 'Cathode';
            itf = 'Interface';

            state.(ct).(itf) = model.updateReactionRateCoefficient(model.(ct).(itf), state.(ct).(itf));
            
        end

        function state = setupElectrolyteMassCons(model, state)

            an    = 'Anode';
            ct    = 'Cathode';
            itf   = 'Interface';
            elyte = 'Electrolyte';

            anArea = model.anodeArea;
            ctArea = model.cathodeArea;

            anR = state.(an).R;
            ctR = state.(ct).(itf).intercalationFlux;
            
            state.(elyte).massCons = anR*anArea + ctR*ctArea;
            
        end
        
        function state = updateControl(model, state, drivingForces)
        % note : The variables ctrlVal and ctrlType, which are updated here, are not AD variables
            ctrl = "Control";
            
            switch model.(ctrl).controlPolicy
              case {'CCCV', 'CV'}

                % nothing to do here
                
              case 'CCDischarge'
                
                E    = state.(ctrl).E;
                I    = state.(ctrl).I;
                time = state.time;
                
                [ctrlVal, ctrltype] = drivingForces.src(time, value(I), value(E));
                
                state.(ctrl).ctrlVal  = ctrlVal;
                state.(ctrl).ctrlType = ctrltype;
              otherwise
                error('control policy not recognized');
            end
            
        end



        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            
            ctrl = 'Control';
            switch model.(ctrl).controlPolicy
              case 'CCCV'
                forces.CCCV = true;
              case 'CCDischarge'
                forces.CCDischarge = true;
                forces.src = [];
              case 'CV'
                forces.CV = true;
                forces.src = [];
              otherwise
                error('Error controlPolicy not recognized');
            end
            %NB hack to get thing go
            forces.Imax = [];%model.Control.Imax;
            
        end


        function [state, report] = updateState(model, state, problem, dx, forces)
        
            [state, report] = updateState@BaseModel(model, state, problem, dx, forces);
            
            ctrl = 'Control';            
            state.(ctrl) = model.(ctrl).updateControlState(state.(ctrl));
            
        end
        
        
        function model = validateModel(model, varargin)
        % we do nothing here
        end
        
        function [convergence, values, names] = checkConvergence(model, problem, varargin)
            
            [convergence, values, names] = checkConvergence@BaseModel(model, problem, varargin{:});
            if all(convergence)

                ctrl = 'Control';            
                
                Emin = model.(ctrl).lowerCutoffVoltage;
                Imax = model.(ctrl).Imax;
                
                state = problem.state;

                ctrlType = state.(ctrl).ctrlType;
                E = value(state.(ctrl).E);
                I = value(state.(ctrl).I);
                
                if (strcmp(ctrlType, 'CC_discharge1') & (E <= Emin)) || (strcmp(ctrlType, 'CV_discharge2') & (I > Imax))
                    convergence(:) = Inf;
                end
            
            end
            
        end
        
        
        function cleanState = addStaticVariables(model, cleanState, state, state0)

            cleanState = addStaticVariables@BaseModel(model, cleanState, state);

            an   = 'Anode';
            ct   = 'Cathode';
            sei  = 'SolidElectrodeInterface';
            itf  = 'Interface';
            sr   = 'SideReaction';
            ctrl = 'Control';
            
            cleanState.T = state.T;
            % zero potential boundary condition at anode
            cleanState.(an).phi = 0;
            % zero external potential drop at the cathode
            cleanState.(ct).(itf).externalPotentialDrop = 0;
            % external EC concentration is set to constant
            cleanState.(an).(sei).cExternal = state.(an).(sei).cExternal;
            
            cleanState.(ctrl) = model.(ctrl).addStaticVariables(cleanState.(ctrl), state.(ctrl));
                
        end
        
        
        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            
            [model, state] = prepareTimestep@BaseModel(model, state, state0, dt, drivingForces);
            
            ctrl = 'Control';
            state.(ctrl) = model.(ctrl).prepareStepControl(state.(ctrl), state0.(ctrl), dt, drivingForces);
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            
            [state, report] = updateAfterConvergence@BaseModel(model, state0, state, dt, drivingForces);
            
            ctrl = 'Control';
            state.(ctrl) = model.(ctrl).updateControlAfterConvergence(state.(ctrl), state0.(ctrl), dt);
        end
        
    end

    methods(Static)
        
        function state = updateReactionRateCoefficient(interfaceModel, state)

            cmax = interfaceModel.cmax;
            k0   = interfaceModel.k0;
            n    = interfaceModel.n;
            R    = interfaceModel.constants.R;
            F    = interfaceModel.constants.F;

            c      = state.cElectrodeSurface;
            
            % We use regularizedSqrt to regularize the square root function and avoid the blow-up of derivative at zero.
            th = 1e-3*cmax;
            j0 = k0.*regularizedSqrt((cmax - c).*c, th)*n*F;

            state.j0 = j0;


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
