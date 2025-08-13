classdef Interface < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();


        %% Input parameters

        % Standard parameters

        saturationConcentration      % the saturation concentration of the guest molecule in the host material
        numberOfElectronsTransferred % stoichiometric number of electrons transferred in the electrochemical reaction
        volumetricSurfaceArea        % surface area of the active material - electrolyte interface per volume of electrode
        activationEnergyOfReaction   % the activation energy of the electrochemical reaction
        reactionRateConstant         % the reaction rate constant of the electrochemical reaction

        % The exchange current at the electrode-electrolyte interface under equilibrium conditions.
        % Tf empty, the default expression using the reaction rate constant is used, see method
        % Interface.updateReactionRateCoefficient. The function is given as a struct with the fields:
        %   - type = {"function", "constant"} % if "constant" is selected, we use the reactionRateConstant value
        %   - functionName :  matlab function name (should be available in path)
        %   - argumentList = ["cElectrodeSurface", "cmax"]
        exchangeCurrentDensity

        guestStoichiometry100 % the ratio of the concentration of the guest molecule to the saturation concentration
                              % of the guest molecule in a phase at a cell voltage that is defined as 100% SOC
        guestStoichiometry0   % the ratio of the concentration of the guest molecule to the saturation concentration
                              % of the guest molecule in a phase at a cell voltage that is defined as 0% SOC
        density               % the mass density of the active material

        % A function to determine the open-circuit potential of the electrode under given conditions
        % See schema `Utilities/JsonSchemas/Function.schema.json` for a complete description of the function interface
        openCircuitPotential

        % A function to determine the entropy change
        % See schema `Utilities/JsonSchemas/Function.schema.json` for a complete description of the function interface        
        entropyChange

        referenceTemperature % Used to compute effective OCP from reference OCP and entropy change
        
        chargeTransferCoefficient % the charge transfer coefficient that enters in the Butler-Volmer equation (symbol: alpha)

        %% Double layer capacity
        % We use modeling equation from
        % @article{Legrand_2014, title={Including double-layer capacitance in lithium-ion battery mathematical models},
        % journal={Journal of Power Sources},
        % author={Legrand, N. and Raël, S. and Knosp, B. and Hinaje, M. and Desprez, P. and Lapicque, F.}, year={2014}}

        useDoubleLayerCapacity % if true, add double layer capacity (default is false)
        doubleLayerCapacitance % Value of electric double layer capacitance / Fm^-2

        %% Computed parameters at model setup

        computeOCPFunc           % Function object to compute the OCP (see Utilities/FunctionInterface/Function.m)
        computeOCP               % function handler to compute the OCP (can be called directly)
        includeEntropyChange     % flag which determines if entropy change should be computed and included
        computeEntropyChangeFunc % Function object to compute dUdT (see Utilities/FunctionInterface/Function.m)
        computeEntropyChange     % function handler to compute dUdT (can be called directly)
        useJ0Func                % true if we use a function to compute the function computeJ0Func to compute the exchange current density
        computeJ0Func            % used when useJ0Func is true. Function object to compute J0 as function of cElectrode, see method updateReactionRateCoefficient
        computeJ0                % function handler used when useJ0Func is true. (can be called directly)

    end

    methods

        function model = Interface(inputparams)

            model = model@BaseModel();

            fdnames = {'G'                           , ...
                       'saturationConcentration'     , ...
                       'numberOfElectronsTransferred', ...
                       'volumetricSurfaceArea'       , ...
                       'activationEnergyOfReaction'  , ...
                       'reactionRateConstant'        , ...
                       'exchangeCurrentDensity'      , ...
                       'guestStoichiometry100'       , ...
                       'guestStoichiometry0'         , ...
                       'density'                     , ...
                       'openCircuitPotential'        , ...
                       'entropyChange'               , ...
                       'referenceTemperature'        , ...
                       'includeEntropyChange'        , ...
                       'chargeTransferCoefficient'   , ...
                       'useDoubleLayerCapacity'      , ...
                       'doubleLayerCapacitance'};

            model = dispatchParams(model, inputparams, fdnames);

            [model.computeOCP, ...
             model.computeOCPFunc] = setupFunction(inputparams.openCircuitPotential);

            if model.includeEntropyChange
                [model.computeEntropyChange, ...
                 model.computeEntropyChangeFunc] = setupFunction(inputparams.entropyChange);
            end
            
            j0 = inputparams.exchangeCurrentDensity;

            if ~isempty(j0)
                if isa(j0, 'struct') && isfield(j0, 'functionFormat')
                    model.useJ0Func = true;
                    [model.computeJ0, ...
                     model.computeJ0Func] = setupFunction(inputparams.exchangeCurrentDensity);
                else
                    model.useJ0Func = false;
                end
            else
                model.useJ0Func = false;
            end

            if isempty(model.useDoubleLayerCapacity)
                model.useDoubleLayerCapacity = false;
            end
                
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Temperature
            varnames{end + 1} = 'T';
            % potential in electrode
            varnames{end + 1} = 'phiElectrode';
            % charge carrier concentration in electrode - value at surface
            varnames{end + 1} = 'cElectrodeSurface';
            % potential in electrolyte
            varnames{end + 1} = 'phiElectrolyte';
            % charge carrier concentration in electrolyte
            varnames{end + 1} = 'cElectrolyte';
            % Electrode over potential
            varnames{end + 1} = 'eta';
            % Reaction rate [mol s^-1 m^-2]
            varnames{end + 1} = 'R';
            % OCP0 [V] at reference temperature
            varnames{end + 1} = 'OCP0';
            % OCP [V]
            varnames{end + 1} = 'OCP';
            % Entropy
            if model.includeEntropyChange
                varnames{end + 1} = 'dUdT';
            end
            % Reaction rate coefficient [A m^-2]
            varnames{end + 1} = 'j0';

            model = model.registerVarNames(varnames);

            if model.useDoubleLayerCapacity
                varnames = {};
                % Double layer capacity rate [mol s^-1 m^-2]
                varnames{end + 1} = 'capacityR';
                % Reaction rate [mol s^-1 m^-2] (same as R without the double layer capacity, now R will contain the sum)
                varnames{end + 1} = 'reactionR';
                % Double layer capacity rate equation
                varnames{end + 1} = 'capacityRequation';
                model = model.registerVarNames(varnames);                
            end
            
            fn = @Interface.updateReactionRateCoefficient;
            
            if model.useJ0Func
                switch computeJ0.numberOfArguments
                  case 1                
                    inputnames = {'cElectrodeSurface'};
                  case 3
                    inputnames = {'cElectrodeSurface', 'T', 'cElectrolyte'};
                  otherwise
                    error('number of argument not supported')
                end
            else
                inputnames = {'T', 'cElectrolyte', 'cElectrodeSurface'};
            end
            model = model.registerPropFunction({'j0', fn, inputnames});

            fn = @Interface.updateOCP0;
            inputnames = {'cElectrodeSurface'};
            model = model.registerPropFunction({'OCP0', fn, inputnames});

            if model.includeEntropyChange
                fn = @Interface.updatedUdT;
                inputnames = {'cElectrodeSurface'};            
                model = model.registerPropFunction({'dUdT', fn, inputnames});
            end

            fn = @Interface.updateOCP;
            if model.includeEntropyChange
                inputnames = {'OCP0', 'dUdT', 'T'};
            else
                inputnames = {'OCP0'};
            end
            model = model.registerPropFunction({'OCP', fn, inputnames});
            
            fn = @Interface.updateEta;
            inputnames = {'phiElectrolyte', 'phiElectrode', 'OCP'};
            model = model.registerPropFunction({'eta', fn, inputnames});

            % This function is used when SEI layer
            % fn = @Interface.updateEtaWithEx;
            % inputnames = {'phiElectrolyte', 'phiElectrode', 'OCP', 'externalPotentialDrop'};
            % model = model.registerPropFunction({'eta', fn, inputnames});


            if ~model.useDoubleLayerCapacity            

                fn = @Interface.updateReactionRate;
                inputnames = {'T', 'eta', 'j0'};
                model = model.registerPropFunction({'R', fn, inputnames});
                
            else
                
                fn = @Interface.updateReactionCapacityRateEquation;
                fn = {fn, @(prop) PropFunction.accumFuncCallSetupFn(prop)};
                inputnames = {'phiElectrolyte', 'phiElectrode', 'cElectrolyte', 'capacityR', 'T'};
                model = model.registerPropFunction({'capacityRequation', fn, inputnames});

                fn = @Interface.updateReactionRateWithCapacity;
                inputnames = {'T', 'eta', 'j0'};
                model = model.registerPropFunction({'reactionR', fn, inputnames});
                
                fn = @Interface.updateTotalRateWithCapacity;
                inputnames = {'reactionR', 'capacityR'};
                model = model.registerPropFunction({'R', fn, inputnames});
                
            end


        end

        function jsonstruct = exportParams(model)

            jsonstruct = exportParams@BaseModel(model);

            fdnames = {'saturationConcentration'     , ...
                       'numberOfElectronsTransferred', ...
                       'volumetricSurfaceArea'       , ...
                       'activationEnergyOfReaction'  , ...
                       'reactionRateConstant'        , ...
                       'exchangeCurrentDensity'      , ...
                       'guestStoichiometry100'       , ...
                       'guestStoichiometry0'         , ...
                       'openCircuitPotential'        , ...
                       'chargeTransferCoefficient'};

            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                jsonstruct.(fdname) = model.(fdname);
            end


        end

        function state = dipatchTemperature(model, state)

            sd = 'SolidDiffusion';
            state.(sd).T = state.T;

        end

        function state = updateOCP0(model, state)
            
            computeOCPFunc = model.computeOCPFunc;
            cmax       = model.saturationConcentration;
            
            c = state.cElectrodeSurface;

            switch computeOCPFunc.numberOfArguments
              case 1
                state.OCP0 = computeOCPFunc.eval(c/cmax);
              otherwise
                error('number of argument not recognized')
            end

        end

        function state = updatedUdT(model, state)

            computeEntCh = model.computeEntropyChangeFunc;
            cmax         = model.saturationConcentration;

            c = state.cElectrodeSurface;
            
            if computeEntCh.numberOfArguments
                state.dUdT = computeEntCh(c/cmax);
            else
                state.dUdT = computeEntCh(c, cmax);
            end

        end
        
        function state = updateOCP(model, state)

            if model.includeEntropyChange
                
                Tref = model.referenceTemperature;
                
                dUdT = state.dUdT;
                T    = state.T;
                OCP0 = state.OCP0;

                state.OCP = state.OCP0 + dUdT.*(T - Tref);
                
            else
                
                state.OCP = state.OCP0;
                
            end
            
        end

        function state = updateReactionRateCoefficient(model, state)

            if model.useJ0Func

                computeJ0 = model.computeJ0Func;

                switch computeJ0.numberOfArguments

                    case 1

                    cmax      = model.saturationConcentration;
                    theta0    = model.guestStoichiometry0;
                    theta100  = model.guestStoichiometry100;
                    
                    c = state.cElectrodeSurface;
                    
                    cmin = theta0*cmax;
                    cmax = theta100*cmax;
                    
                    soc = (c - cmin)./(cmax - cmin);
                    
                    j0 = computeJ0(soc);
                    
                  case 4

                    cmax = model.saturationConcentration;
                    
                    celyte = state.cElectrolyte;
                    celde  = state.cElectrodeSurface;
                    T      = state.T;

                    j0 = computeJ0(celyte, celde, cmax, T);

                  otherwise

                    error('number of argument not recgonized');
                    
                end

            else

                Tref = 298.15;  % [K]

                cmax = model.saturationConcentration;
                k0   = model.reactionRateConstant;
                Eak  = model.activationEnergyOfReaction;
                n    = model.numberOfElectronsTransferred;
                F    = model.constants.F;
                R    = model.constants.R;

                T      = state.T;
                cElyte = state.cElectrolyte;
                c      = state.cElectrodeSurface;

                % Calculate reaction rate constant
                k = k0.*exp(-Eak./R.*(1./T - 1/Tref));

                % We use regularizedSqrt to regularize the square root function and avoid the blow-up of derivative at zero.
                th = 1e-3*cmax;
                coef = cElyte.*(cmax - c).*c;
                coef(coef < 0) = 0;
                j0 = k.*regularizedSqrt(coef, th)*n*F;

            end

            state.j0 = j0;

        end

        
        function state = updateEta(model, state)

            phiElyte = state.phiElectrolyte;
            phiElde  = state.phiElectrode;
            OCP      = state.OCP;

            state.eta = (phiElde - phiElyte - OCP);

        end

        function state = updateEtaWithEx(model, state)

            phiElyte = state.phiElectrolyte;
            phiElde  = state.phiElectrode;
            OCP      = state.OCP;
            dphi     = state.externalPotentialDrop;

            state.eta = (phiElde - phiElyte - OCP - dphi);

        end

        function R = computeRate(model, state)
        % From definition of the overpotential eta, we have that reaction rate R is positive for oxydation.

            n     = model.numberOfElectronsTransferred;
            F     = model.constants.F;
            alpha = model.chargeTransferCoefficient;

            T   = state.T;
            j0  = state.j0;
            eta = state.eta;

            R = ButlerVolmerEquation(j0, alpha, n, eta, T);

            R = R/(n*F); % reaction rate in mol/(s*m^2)
            
        end

        function state = updateReactionRate(model, state)

            state.R = model.computeRate(state);
            
        end
        
        function state = updateReactionCapacityRateEquation(model, state, state0, dt)

            cDL = model.doubleLayerCapacitance;
            F   = model.constants.F;
            R   = model.constants.R;
            
            jDL   = state.capacityR;
            T     = state.T;
            c     = state.cElectrolyte;
            c0    = state0.cElectrolyte;
            dphi  = state.phiElectrode - state.phiElectrolyte;
            dphi0 = state0.phiElectrode - state0.phiElectrolyte;

            state.capacityRequation = jDL - (cDL/(F*dt))*((dphi - dphi0) + (R.*T./F)./c.*(c - c0));
            
        end
        
        function state = updateReactionRateWithCapacity(model, state)

            state.reactionR = model.computeRate(state);
            
        end
        
        function state = updateTotalRateWithCapacity(model, state)

            state.R = state.reactionR + state.capacityR;
        end

        function newstate = addVariablesAfterConvergence(model, newstate, state)

            if model.useDoubleLayerCapacity
                
                newstate.cElectrolyte   = state.cElectrolyte;
                newstate.phiElectrode   = state.phiElectrode;
                newstate.phiElectrolyte = state.phiElectrolyte;
                
            end
        
        end
        
    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes



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
