classdef ActiveMaterial < BaseModel
    
    properties
        
        %
        % instance of :class:`Interface <Electrochemistry.Electrodes.Interface>`
        %

        %% Sub-Models
        
        Interface
        SolidDiffusion        

        %% Input parameters

        % Standard parameters

        electronicConductivity % the electronic conductivity of the material (symbol: sigma)
        density                % the mass density of the material (symbol: rho)
        massFraction           % the ratio of the mass of the material to the total mass of the phase or mixture (symbol: gamma)
        
        thermalConductivity    % the intrinsic Thermal conductivity of the active component
        specificHeatCapacity   % Specific Heat capacity of the active component

        diffusionModelType     % either 'full' or 'simple'
        
        % Advanded parameters
        
        standAlone % Set to true if model is used as main model for development purpose (standard use is as a sub-model)

        % Coupling parameters
        
        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)


        % Others
        
        % properties needed when model is run in stand-alone mode
        
        funcCallList
        primaryVarNames
        equationVarNames
        equationNames
        equationTypes
        
    end
    
    methods
        
        function model = ActiveMaterial(paramobj)
        %
        % ``paramobj`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %
            model = model@BaseModel();
            
            fdnames = {'electronicConductivity', ... 
                       'density'               , ...                
                       'massFraction'          , ...
                       'thermalConductivity'   , ...    
                       'specificHeatCapacity'  , ...   
                       'volumeFraction'        , ... 
                       'externalCouplingTerm'  , ...
                       'diffusionModelType'    , ...
                       'standAlone'};

            model = dispatchParams(model, paramobj, fdnames);

            model.Interface = Interface(paramobj.Interface);

            diffusionModelType = model.diffusionModelType;

            switch model.diffusionModelType
              case 'simple'
                model.SolidDiffusion = SimplifiedSolidDiffusionModel(paramobj.SolidDiffusion);
              case 'full'
                model.SolidDiffusion = FullSolidDiffusionModel(paramobj.SolidDiffusion);
              otherwise
                error('Unknown diffusionModelType %s', diffusionModelType);
            end

            if model.standAlone

                model = model.setupStandAloneModel();
                
            end
            
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            varnames = {'T'};
            model = model.registerVarNames(varnames);

            fn = @ActiveMaterial.dispatchTemperature;
            model = model.registerPropFunction({{sd, 'T'}, fn, {'T'}});
            model = model.registerPropFunction({{itf, 'T'}, fn, {'T'}});
            
            if model.standAlone

                varnames = {};

                % Volumetric current in A/m^3. It corresponds to the current density multiplied with the volumetric
                % surface area.
                varnames{end + 1} = 'I';
                % Potential at Electrode
                varnames{end + 1} = 'E';
                % Charge Conservation equation
                varnames{end + 1} = 'chargeCons';
                
                model = model.registerVarNames(varnames);

                varnames = {{itf, 'dUdT'}, ...
                            'jCoupling', ...
                            'jExternal'};
                model = model.removeVarNames(varnames);

                varnames = {'T', ...
                            {itf, 'cElectrolyte'},... 
                            {itf, 'phiElectrolyte'}};
                model = model.registerStaticVarNames(varnames);
                
            end

            fn = @ActiveMaterial.updateRvol;
            model = model.registerPropFunction({{sd, 'Rvol'}, fn, {{itf, 'R'}}});
            
            fn = @ActiveMaterial.updateConcentrations;
            model = model.registerPropFunction({{itf, 'cElectrodeSurface'}, fn, {{sd, 'cSurface'}}});
            
            if model.standAlone
                
                fn = @ActiveMaterial.updateControl;
                fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
                model = model.registerPropFunction({'I', fn, {}});
                
                fn = @ActiveMaterial.updateChargeCons;
                inputnames = {'I', ...
                              {sd, 'Rvol'}};
                model = model.registerPropFunction({'chargeCons', fn, inputnames});

                fn = @ActiveMaterial.updatePhi;
                model = model.registerPropFunction({{itf, 'phiElectrode'}, fn, {'E'}});

            end
            
        end
        

        function model = setupStandAloneModel(model)

            model = model.setupComputationalGraph();

            cgt = model.computationalGraph();
            
            model.funcCallList     = cgt.getOrderedFunctionCallList();
            model.primaryVarNames  = cgt.getPrimaryVariableNames();
            model.equationVarNames = cgt.getEquationVariableNames();
            
            function str = setupName(varname)
                shortvarname = cellfun(@(elt) Battery.shortenName(elt), varname, 'uniformoutput', false);
                str = Battery.varToStr(shortvarname);
            end
            model.equationNames = cellfun(@(varname) setupName(varname), model.equationVarNames, 'uniformoutput', false);
            model.equationTypes = repmat({'cell'}, 1, numel(model.equationNames));
            
        end

        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            time = state0.time + dt;
            state = model.initStateAD(state);
            
            %% We call the assembly equations ordered from the graph

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            % some scaling of the equations

            %% Setup equations and add some scaling
            n  = model.(itf).numberOfElectronsTransferred; % number of electron transfer (equal to 1 for Lithium)
            F  = model.(sd).constants.F;
            rp = model.(sd).particleRadius;
            
            scalingcoef = n*F/(4*pi*rp^3/3);

            state.(sd).massCons         = scalingcoef*state.(sd).massCons;
            state.(sd).solidDiffusionEq = scalingcoef*state.(sd).solidDiffusionEq;
            
            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end
            
            names       = model.equationNames;
            types       = model.equationTypes;
            primaryVars = model.primaryVarNames;
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end

        function primaryvarnames = getPrimaryVariableNames(model)

            primaryvarnames = model.primaryVarNames;

        end


        function state = updateControl(model, state, drivingForces)
            
            state.I = drivingForces.src(state.time);
            
        end

        function state = updatePhi(model, state)

            itf = 'Interface';
            
            state.(itf).phiElectrode = state.E;
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state, state0)
            
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            
            itf = 'Interface';
            
            cleanState.T = state.T;
            cleanState.I = state.I;
            cleanState.(itf).cElectrolyte   = state.(itf).cElectrolyte;
            cleanState.(itf).phiElectrolyte = state.(itf).phiElectrolyte;
            
        end
        
        function state = updateChargeCons(model, state)
        % Only used for stand-alone model

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            n    = model.(itf).numberOfElectronsTransferred;
            F    = model.(itf).constants.F;
            
            I = state.I;
            Rvol = state.(sd).Rvol;

            state.chargeCons = I - Rvol*n*F;

        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            
            [state, report] = updateAfterConvergence@BaseModel(model, state0, state, dt, drivingForces);
            
        end
         
        function model = validateModel(model, varargin)
        % 
        end

        %% assembly functions use in this model

        function state = updateRvol(model, state)
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            vsa = model.(itf).volumetricSurfaceArea;
            
            Rvol = vsa.*state.(itf).R;
            
            state.(sd).Rvol = Rvol;
            
        end        
        
        function state = updateConcentrations(model, state)

            sd  = 'SolidDiffusion';
            itf = 'Interface';
            
            state.(itf).cElectrodeSurface = state.(sd).cSurface;
            
        end

        
        function state = dispatchTemperature(model, state)

            state.Interface.T      = state.T;
            state.SolidDiffusion.T = state.T;
            
        end


        function state = updateAverageConcentration(model, state)

            % shortcut
            sd  = 'SolidDiffusion';

            vf       = model.volumeFraction;
            am_frac  = model.activeMaterialFraction;
            vols     = model.G.cells.volumes;
            
            c = state.(sd).cAverage;

            vols = am_frac*vf.*vols;

            cAverage = sum(c.*vols)/sum(vols);

            state.cAverage = cAverage;
            
        end
        
        
        
    end

    methods (Static)

        function str = varToStr(varname)

            str = strjoin(varname, '_');

        end

        function str = shortenName(name)

            namemapping = {'ActiveMaterial'   , 'am'  ; ...
                           'ActiveMaterial1'  , 'am1' ; ...
                           'ActiveMaterial2'  , 'am2' ; ...
                           'Interface'        , 'itf' ; ...
                           'SolidDiffusion'   , 'sd'};

            [found, ind] = ismember(name, namemapping(:, 1));

            if found
                str = namemapping{ind, 2};
            else
                str = name;
            end

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
