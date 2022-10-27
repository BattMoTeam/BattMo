classdef Interface < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();

        cmax
        
        % number of electron transfer
        n

        % Physicochemical properties
        volumeFraction
        volumetricSurfaceArea  % Surface area to volume,       [m2 m^-3]
        density                % [kg m^-3]
        theta0                 % Minimum lithiation, 0% SOC    [-]
        theta100               % Maximum lithiation, 100% SOC  [-]
        k0                     % Reference rate constant       [m^2.5 mol^-0.5 s^-1]
        Eak                    % Reaction activation energy    [J mol^-1]
        alpha                  % coefficient in Butler-Volmer equation
        
        computeOCPFunc % Function handler to compute OCP
        
        useJ0Func
        computeJ0Func % used when useJ0Func is true. Function handler to compute J0 as function of cElectrode, see method updateReactionRateCoefficient
        
    end

    methods

        function model = Interface(paramobj)

            model = model@BaseModel();

             % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);

            fdnames = {'G'                    , ...
                       'cmax'                 , ...
                       'n'                    , ...
                       'volumeFraction'       , ...
                       'volumetricSurfaceArea', ...  
                       'density'              , ...                
                       'theta0'               , ...                 
                       'theta100'             , ...               
                       'k0'                   , ...                     
                       'Eak'                  , ...                    
                       'alpha'};

            model = dispatchParams(model, paramobj, fdnames);

            model.computeOCPFunc = str2func(paramobj.OCP.functionname);

            if ~isempty(paramobj.j0)
                switch paramobj.j0.type
                  case 'function'
                    model.useJ0Func = true;
                    model.computeJ0Func = str2func(paramobj.j0.functionname);
                  case 'constant'
                    model.useJ0Func = false;
                  otherwise
                    errror('type of j0 not recognized.')
                end
            else
                model.useJ0Func = false;
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
            % Reaction rate in mol/(s*m^2)
            varnames{end + 1} = 'R';
            % External potential drop used in Butler-Volmer
            % varnames{end + 1} = 'externalPotentialDrop';
            % 
            varnames{end + 1} = 'dUdT';
            % OCP
            varnames{end + 1} = 'OCP';
            % Reaction rate coefficient
            varnames{end + 1} = 'j0';
            
            model = model.registerVarNames(varnames);
            
            fn = @Interface.updateReactionRateCoefficient;
            if model.useJ0Func
                inputnames = {'cElectrodeSurface'};
            else
                inputnames = {'T', 'cElectrolyte', 'cElectrodeSurface'};
            end
            model = model.registerPropFunction({'j0', fn, inputnames});

            fn = @Interface.updateOCP;
            inputnames = {'cElectrodeSurface', 'T'};
            model = model.registerPropFunction({'OCP', fn, inputnames});
            model = model.registerPropFunction({'dUdT', fn, inputnames});
            
            fn = @Interface.updateEta;
            inputnames = {'phiElectrolyte', 'phiElectrode', 'OCP'};            
            model = model.registerPropFunction({'eta', fn, inputnames});
            
            % fn = @Interface.updateEtaWithEx;
            % inputnames = {'phiElectrolyte', 'phiElectrode', 'OCP', 'externalPotentialDrop'};
            % model = model.registerPropFunction({'eta', fn, inputnames});            
            
            fn = @Interface.updateReactionRate;
            inputnames = {'T', 'eta', 'j0'};
            model = model.registerPropFunction({'R', fn, inputnames});
            
            
        end
        
        function state = dipatchTemperature(model, state)

            sd = 'SolidDiffusion';
            state.(sd).T = state.T;
            
        end
        
        function state = updateOCP(model, state)

            computeOCP = model.computeOCPFunc;
            cmax = model.cmax;

            c = state.cElectrodeSurface;
            T = state.T;

            [state.OCP, state.dUdT] = computeOCP(c, T, cmax);
            
        end

        function state = updateReactionRateCoefficient(model, state)


            if model.useJ0Func

                computeJ0 = model.computeJ0Func;
                cmax      = model.cmax;
                theta0    = model.theta0;
                theta100  = model.theta100;
                
                c = state.cElectrodeSurface;

                cmin = theta0*cmax;
                cmax = theta100*cmax;

                soc = (c - cmin)./(cmax - cmin);
                
                j0 = computeJ0(soc);

            else
                
                Tref = 298.15;  % [K]

                cmax = model.cmax;
                k0   = model.k0;
                Eak  = model.Eak;
                n    = model.n;
                F    = model.constants.F;
                R    = model.constants.R;

                T      = state.T;
                cElyte = state.cElectrolyte;
                c      = state.cElectrodeSurface;
                
                % Calculate reaction rate constant
                k = k0.*exp(-Eak./R.*(1./T - 1/Tref));

                % We use regularizedSqrt to regularize the square root function and avoid the blow-up of derivative at zero.
                th = 1e-3*cmax;
                j0 = k.*regularizedSqrt(cElyte.*(cmax - c).*c, th)*n*F;
                
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
        
            
        function state = updateReactionRate(model, state)
        % From definition of the overpotential eta, we have that reaction rate R is positive for oxydation.
            n     = model.n;
            F     = model.constants.F;
            alpha = model.alpha;

            T   = state.T;
            j0  = state.j0;
            eta = state.eta;
            
            R = ButlerVolmerEquation(j0, alpha, n, eta, T);

            state.R = R/(n*F); % reaction rate in mol/(s*m^2)

        end
        
    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes



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
