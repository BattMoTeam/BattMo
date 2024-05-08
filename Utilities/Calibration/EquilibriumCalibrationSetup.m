classdef EquilibriumCalibrationSetup
%% The goal of this class is to solve the following calibration problem:
% Given a discharge curve, find the set of parameters that give the better match.
% By default (calibrationCase = 1), the calibration parameters are the guest stoichiometry at discharge start and the volume fraction, for both electrodes.
    
    properties
        
        model % Battery model

        F % Faraday Constant
        Temperature = 298.15% Temperature (used in OCP curves)

        % Discharge curve measured experimentally at near-equilibrium condition
        exptime   % Time [s]
        expU      % Voltage [V]
        expI      % Current [A]
        totalTime % total time (set as exptime(end))
        
        packingMass % mass of packing 

        bounds % Vector of variable bounds on the variables, with field
               % - lower
               % - upper

        calibrationCase = 1
        % different calibration case depending on the parameters that are chosen. See method printVariableChoice
        % At the moment, the following is implemented
        % case 1 (default) : The calibration parameters are theta100 and volume fraction for both electrodes
        % case 2           : The calibration parameters are theta100 for the negative electrode, the volume fractions for both electrodes
        % case 3           : The calibration parameters are theta100, theta0 and volume fraction for both electrodes and we add a constraint on the np-ratio (thus we use IpOpt solver)
        
        calibrationParameters 

    end

    methods

        function ecs = EquilibriumCalibrationSetup(model, expdata)

            con = PhysicalConstants();
            
            ecs.F         = con.F;
            ecs.exptime   = expdata.time;
            ecs.expU      = expdata.U;
            ecs.expI      = expdata.I;
            ecs.totalTime = ecs.exptime(end);
            ecs.model     = model;

            ecs = ecs.setupCalibrationCase(1, 'verbose', false);
            
        end

        function ecs = setupCalibrationCase(ecs, calibrationCase, varargin)

            opt = struct('verbose', true);
            [opt, extras] = merge_options(opt, varargin{:});

            ecs.calibrationCase = calibrationCase;
            
            switch calibrationCase
                
              case 1

                data_default = struct('ne_theta_max', 1);

              case 2

                data_default = struct('pe_theta', 0.99, ...
                                      'ne_theta_max', 1);
              case 3

                data_default = struct('ne_theta_max', 1, ...
                                      'np_ratio', 1.1);
                
              otherwise
                error('case number not recognized');
            end
            
            ecs.calibrationParameters = merge_options(data_default, extras{:});

            if opt.verbose
                ecs.printVariableChoice();
            end

            ecs = ecs.setupDefaultVariableBounds('verbose', opt.verbose);
            
        end

        function ecs = setupDefaultVariableBounds(ecs, varargin)
            
            opt = struct('verbose', true);
            opt = merge_options(opt, varargin{:});
            
            data = ecs.calibrationParameters;
            
            switch ecs.calibrationCase
                
              case 1
                
                % X(1) : theta100 anode
                % X(2) : volume fraction anode
                % X(3) : theta100 cathode
                % X(4) : volume fraction cathode
                
                % we use a non-zero default value to avoid problem
                bounds.lower = [0.1, 0.1, 0.1, 0.1]';
                bounds.upper = [1, 1, 1, 1]';
                
              case 2
                
                %% recall : ordering of parameters
                %
                % X(1) : theta100 anode
                % X(2) : volume fraction anode
                % X(3) : volume fraction cathode
                
                bounds.lower = [0.1, 0.1, 0.1]';
                bounds.upper = [1, 1, 1]';

              case 3
                
                % X(1) : theta100 anode
                % X(2) : theta0 anode
                % X(3) : volume fraction anode
                % X(4) : theta100 cathode
                % X(5) : theta0 cathode
                % X(6) : volume fraction cathode
                
                % we use a non-zero default value to avoid problem
                bounds.lower = 0.01*ones(6, 1);
                bounds.upper = ones(6, 1);
                                
              otherwise

                error('ecs.calibrationCase not recognized');
                
            end

            if opt.verbose
                fprintf('The variable bounds are\n');
                fprintf('lower\n')
                display(bounds.lower);
                fprintf('upper\n')
                display(bounds.upper);
            end
            
            ecs.bounds = bounds;
            
        end

        function linIneq = setupIneqConstraints(ecs)

            bounds = ecs.bounds;

            n = size(bounds.lower, 1);
            
            linIneq.A = [-eye(n); eye(n)];
            linIneq.b = [bounds.lower; bounds.upper];
            
        end

        function params = setupOptimParams(ecs)
        % For unitBoxBFGS
            linIneq = ecs.setupIneqConstraints();
            params = {'objChangeTol'    , 1e-12, ...
                      'maximize'        , false, ...
                      'maxit'           , 1000 , ...
                      'maxInitialUpdate', 1e-6 , ...
                      'enforceFeasible' , true , ...
                      'lineSearchMaxIt' , 10   , ...
                      'wolfe2'          , 0.99 , ...
                      'linIneq'         , linIneq};
        end
        

        function X = getDefaultValue(ecs)

            
            model = ecs.model;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            switch ecs.calibrationCase

              case 1
                
                %% recall : ordering of parameters
                %
                % X(1) : theta100 anode
                % X(2) : volume fraction anode
                % X(3) : theta100 cathode
                % X(4) : volume fraction cathode
                
                X = nan(4, 1);

                compInds = model.(ne).(co).compInds;
                X(1) = model.(ne).(co).(am).(itf).guestStoichiometry100;
                vf   = model.(ne).(co).volumeFraction*model.(ne).(co).volumeFractions(compInds.(am));
                X(2) = vf;
                
                compInds = model.(pe).(co).compInds;
                X(3) = model.(pe).(co).(am).(itf).guestStoichiometry100;
                vf   = model.(pe).(co).volumeFraction*model.(pe).(co).volumeFractions(compInds.(am));
                X(4) = vf;

              case 2

                %% recall : ordering of parameters
                %
                % X(1) : theta100 anode
                % X(2) : volume fraction anode
                % X(3) : volume fraction cathode
                
                X = nan(3, 1);

                compInds = model.(ne).(co).compInds;
                X(1) = model.(ne).(co).(am).(itf).guestStoichiometry100;
                vf   = model.(ne).(co).volumeFraction*model.(ne).(co).volumeFractions(compInds.(am));
                X(2) = vf;
                
                compInds = model.(pe).(co).compInds;
                vf   = model.(pe).(co).volumeFraction*model.(pe).(co).volumeFractions(compInds.(am));
                X(3) = vf;

              case 3
                
                % X(1) : theta100 anode
                % X(2) : theta0 anode
                % X(3) : volume fraction anode
                % X(4) : theta100 cathode
                % X(5) : theta0 cathode
                % X(6) : volume fraction cathode

                X = nan(6, 1);

                compInds = model.(ne).(co).compInds;
                X(1) = model.(ne).(co).(am).(itf).guestStoichiometry100;
                X(2) = model.(ne).(co).(am).(itf).guestStoichiometry0;
                vf   = model.(ne).(co).volumeFraction*model.(ne).(co).volumeFractions(compInds.(am));
                X(3) = vf;
                
                compInds = model.(pe).(co).compInds;
                X(4) = model.(pe).(co).(am).(itf).guestStoichiometry100;
                X(5) = model.(pe).(co).(am).(itf).guestStoichiometry0;
                vf   = model.(pe).(co).volumeFraction*model.(pe).(co).volumeFractions(compInds.(am));
                X(6) = vf;
                
              otherwise
                error('calibrationCase not recognized')
            end
        end

        function printVariableChoice(ecs)

            switch ecs.calibrationCase
              case 1
                fprintf('\nThe calibration parameters are theta100 and volume fraction for both electrodes\n');
              case 2
                fprintf('\nThe calibration parameters are theta100 for the negative electrode\n');
                fprintf('The volume fractions for both electrodes\n');
              case 3
                fprintf('\nThe calibration parameters are theta100, theta0, and volume fraction for both electrodes\n');
                fprintf('In addition we have a given np_ratio as a constraint\n');
              otherwise
                error('calibrationCase not recognized');
            end
            
               
        end
        
        function vals = getPhysicalValues(ecs, X)

            model = ecs.model;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};
            
            switch ecs.calibrationCase

              case 1
                 
                %% recall : ordering of parameters
                %
                % X(1) : theta100 anode
                % X(2) : volume fraction anode
                % X(3) : theta100 cathode
                % X(4) : volume fraction cathode

                for ielde = 1 : numel(eldes)

                    elde = eldes{ielde};
                    
                    vol  = sum(model.(elde).(co).G.getVolumes());
                    cmax = model.(elde).(co).(am).(itf).saturationConcentration;

                    vals.(elde).theta100       = X(2*ielde - 1);
                    vals.(elde).alpha          = X(2*ielde)*vol*cmax;
                    vals.(elde).volumeFraction = X(2*ielde);
                    
                end

              case 2
                
                %% recall : ordering of parameters
                %
                % X(1) : theta100 anode
                % X(2) : volume fraction anode
                % X(3) : volume fraction cathode


                vol  = sum(model.(ne).(co).G.getVolumes());
                cmax = model.(ne).(co).(am).(itf).saturationConcentration;
                
                vals.(ne).theta100       = X(1);
                vals.(ne).volumeFraction = X(2);
                vals.(ne).alpha          = X(2)*vol*cmax;

                vol  = sum(model.(pe).(co).G.getVolumes());
                cmax = model.(pe).(co).(am).(itf).saturationConcentration;
                
                vals.(pe).volumeFraction = X(3);
                vals.(pe).alpha          = X(3)*vol*cmax;

              case 3
                                
                % X(1) : theta100 anode
                % X(2) : theta0 anode
                % X(3) : volume fraction anode
                % X(4) : theta100 cathode
                % X(5) : theta0 cathode
                % X(6) : volume fraction cathode


                for ielde = 1 : numel(eldes)

                    elde = eldes{ielde};
                    
                    vol  = sum(model.(elde).(co).G.getVolumes());
                    cmax = model.(elde).(co).(am).(itf).saturationConcentration;

                    vals.(elde).theta100       = X(3*ielde - 2);
                    vals.(elde).theta0         = X(3*ielde - 1);
                    vals.(elde).alpha          = X(3*ielde)*vol*cmax;
                    vals.(elde).volumeFraction = X(3*ielde);
                    
                end


              otherwise

                error('calibration case not recognized');
                
            end
                
        end
        
        
        function [fexp, fcomp] = setupfunction(ecs)

            fcomp = @(t, X) ecs.computeF(t, X);
            fexp  = @(t) ecs.experimentalF(t);
            
        end


        function theta = conc(ecs, t, elde, tf, theta_tf, alpha)
        % theta_tf is lithiation at time tf*totalTime
        % returns lithiation at time t
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';

            I = ecs.expI;
            F = ecs.F;
            T = ecs.totalTime;
            
            switch elde
              case ne
                sgn = -1;
              case pe
                sgn = 1;
            end

            theta = theta_tf + (t - tf*T)*((sgn*I)./(F*alpha));
            
        end

        function f = computeF(ecs, t, X)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            thermal = 'ThermalModel';
            ctrl    = 'Control';
            
            I     = ecs.expI;
            model = ecs.model;
            T     = ecs.Temperature;

            vals = ecs.computeThetas(X);

            theta100 = vals.(pe).theta100;
            alpha    = vals.(pe).alpha;
            
            theta = ecs.conc(t, pe, 0, theta100, alpha);
            f = model.(pe).(co).(am).(itf).computeOCPFunc(theta, T, 1);

            theta100 = vals.(ne).theta100;
            alpha    = vals.(ne).alpha;
            
            theta = ecs.conc(t, ne, 0, theta100, alpha);
            f = f - model.(ne).(co).(am).(itf).computeOCPFunc(theta, T, 1);

        end

        function f = experimentalF(ecs, t);
            
            time = ecs.exptime; % in second
            U    = ecs.expU;

            f = interp1(time, U, t, 'linear', 'extrap');
            
        end

        function [z, dz] = objective(ecs, X, varargin)
            
            opt = struct('t', []);
            opt = merge_options(opt, varargin{:});

            t  = opt.t(:); % column vector

            [fexp, fcomp] = ecs.setupfunction();
            X = initVariablesADI(X);

            fdiff = (fcomp(t, X) - fexp(t)).^2;

            g = 0.5 * sum(diff(t) .* (fdiff(1:end-1) + fdiff(2:end)));

            z  = g.val;
            dz = g.jac{1}';
            
        end

        function z = objective_func(ecs, X)
        % for ipopt api
            
            z = ecs.objective(X);
            
        end
        
        function dz = gradient_func(ecs, X)
        % for ipopt api
            
            [~, dz] = ecs.objective(X);
            
        end

        
        function np_ratio = computeNPratio(ecs, X)

            vals = ecs.computeThetas(X);
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            
            np_ratio = (vals.(ne).alpha*(vals.(ne).theta0 - vals.(ne).theta100))./(vals.(pe).alpha*(vals.(pe).theta100 - vals.(pe).theta0));
            
        end
        
        function y = constraints_func(ecs, X)
        % for ipopt api

            assert(ecs.calibrationCase == 3, sprintf('this function has not been implemented for calibration case %d', ecs.calibrationCase));

            data = ecs.calibrationParameters;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            np_ratio = ecs.computeNPratio(X);
            
            y = np_ratio./data.np_ratio;
            
        end

        function dy = jacobian_constraints_func(ecs, X)
        % for ipopt api
            
            X = initVariablesADI(X);
            y = ecs.constraints_func(X);
            dy = y.jac{1};
            
        end

        function y = jacobian_constraints_structure_func(ecs)
        % we do not bother here and setup as full
            
            assert(ecs.calibrationCase == 3, sprintf('this function has not been implemented for calibration case %d', ecs.calibrationCase));
            y = sparse(ones(1, 6));
            
        end

        function vals = computeThetas(ecs, X)

            totalTime = ecs.totalTime;
            vals = ecs.getPhysicalValues(X);

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            eldes = {ne, pe};
            
            switch ecs.calibrationCase

              case 1

                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    vals.(elde).theta0 = ecs.conc(totalTime, elde, 0, vals.(elde).theta100, vals.(elde).alpha);
                end

              case 2

                data = ecs.calibrationParameters;
                vals.(pe).theta100 = data.pe_theta;

                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    vals.(elde).theta0 = ecs.conc(totalTime, elde, 0, vals.(elde).theta100, vals.(elde).alpha);
                end
                
              case 3
                
                % nothing to do
                
              otherwise
                
                error('calibration case not recognized');
                
            end

        end

        
        function vals = computeCapacities(ecs, X)
            
            model = ecs.model;
            F         = ecs.F;
            T         = ecs.Temperature;
            totalTime = ecs.totalTime;
            
            vals = ecs.computeThetas(X);
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};
                
                cmax  = model.(elde).(co).(am).(itf).saturationConcentration;
                
                theta0   = vals.(elde).theta0;
                theta100 = vals.(elde).theta100;
                alpha    = vals.(elde).alpha;
                
                switch elde
                  case ne
                    cM = theta100*cmax;
                    cm = theta0*cmax;
                  case pe
                    cM = theta0*cmax;
                    cm = theta100*cmax;
                end

                cap = (cM - cm)/cmax*alpha*F;
                
                vals.(elde).cmax     = cmax;
                vals.(elde).cap      = cap;
                
            end

        end
            
        
        function props = computeProperties(ecs, X)

            model = ecs.model;
            T     = ecs.Temperature;
            
            props = ecs.computeCapacities(X);
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            
            props.cap = min(props.(ne).cap, props.(pe).cap);
            
            props.(ne).r = props.cap/props.(ne).cap;
            props.(pe).r = props.cap/props.(pe).cap;

            N = 1000; % discretization parameter

            eldes = {ne, pe};
            
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                
                smax = props.(elde).r;
                cmax = props.(elde).cmax;
                c0   = props.(elde).theta100*cmax;
                cT   = props.(elde).theta0*cmax;
                cap  = props.(elde).cap;

                s = smax.*linspace(0, 1, N + 1)';
                c = (1 - s).*c0 + s.*cT;
                f = model.(elde).(co).(am).(itf).computeOCPFunc(c(1 : end - 1), 298, cmax);

                props.(elde).dischargeFunc = @(s) model.(elde).(co).(am).(itf).computeOCPFunc((1 - s).*c0 + s.*cT, T, cmax);
                
                energies.(elde) = cap*smax/N*sum(props.(elde).dischargeFunc(s));

                props.(elde).U = model.(elde).(co).(am).(itf).computeOCPFunc(c0, T, cmax);
                
            end

            props.energy = energies.(pe) - energies.(ne);

            mass = ecs.getMass();
            props.specEnergy = props.energy/mass;

            % props.NPratio = props.(ne).cap/props.(pe).cap;

            props.U = props.(pe).U - props.(ne).U;

        end

        function mass = getMass(ecs)
            
            model       = ecs.model;
            packingMass = ecs.packingMass;
            if isempty(packingMass)
                fprintf('Packing Mass is not given, we set it equal to zero\n');
                packingMass = 0;
            end
            mass = computeCellMass(model, 'packingMass', packingMass);

        end


        function [Xopt, hist] = runUnitBoxBFGS(ecs)

            [fexp, fcomp] = ecs.setupfunction();

            X0 = ecs.getDefaultValue();
            exptime = ecs.exptime;

            topt = exptime;

            opts = {'t', topt};
            f = @(X) ecs.objective(X, opts{:});

            params = ecs.setupOptimParams();

            [~, Xopt, hist] = unitBoxBFGS(X0, f, params{:});

        end

        function [Xopt, info] = runIpOpt(ecs, ipopt_options)
            
            [fexp, fcomp] = ecs.setupfunction();

            X0 = ecs.getDefaultValue();

            funcs.objective         = @(x) ecs.objective_func(x);
            funcs.gradient          = @(x) ecs.gradient_func(x);
            funcs.constraints       = @(x) ecs.constraints_func(x);            
            funcs.jacobian          = @(x) ecs.jacobian_constraints_func(x);
            funcs.jacobianstructure = @()  ecs.jacobian_constraints_structure_func();
            
            options.lb = ecs.bounds.lower;
            options.ub = ecs.bounds.upper;

            assert(ecs.calibrationCase == 3, 'the constraint should be only set ip for calibration case 3');
            
            options.cl = 0;
            options.cu = 0;

            if nargin > 1
                options.ipopt = ipopt_options;
            end
            options.ipopt.hessian_approximation = 'limited-memory';
            
            [Xopt, info] = ipopt(X0, funcs, options);
            
        end

        function printParameters(ecs, X)

            vals = ecs.getPhysicalValues(X);
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';

            eldes = {ne, pe};

            switch ecs.calibrationCase
              case 1
                fprintf('%-25s%20s%20s\n', '', 'theta100', 'volume fraction');
                thetastr = sprintf('%6.5f', vals.(ne).theta100);
                fprintf('%-25s%20s%20.5f \n', ne, thetastr, vals.(ne).volumeFraction);
                thetastr = sprintf('%6.5f', vals.(pe).theta100);
                fprintf('%-25s%20s%20.5f \n', pe, thetastr, vals.(pe).volumeFraction);
              case 2
                fprintf('%-25s%20s%20s\n', '', 'theta100', 'volume fraction');
                thetastr = sprintf('%6.5f', vals.(ne).theta100);
                fprintf('%-25s%20s%20.5f \n', ne, thetastr, vals.(ne).volumeFraction);
                fprintf('%-25s%20s%20.5f \n', pe, '', vals.(pe).volumeFraction);
              case 3
                fprintf('%-25s%20s%20s%20s\n', '', 'theta100', 'theta0', 'volume fraction');
                fprintf('%-25s%20.5f%20.5f%20.5f \n', ne, vals.(ne).theta100, vals.(ne).theta0, vals.(ne).volumeFraction);
                fprintf('%-25s%20.5f%20.5f%20.5f \n', pe, vals.(pe).theta100, vals.(pe).theta0, vals.(pe).volumeFraction);
              otherwise
                error('calibration case not recognized');
            end
            
        end

        function printProperties(ecs, X)

            props = ecs.computeProperties(X);

            mass = ecs.getMass();
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';

            fprintf('%20s: %g [V]\n','initial voltage', props.U);
            fprintf('%20s: %g [Ah]\n', 'Capacity', props.cap/(1*hour));
            fprintf('%20s: %g [Wh/kg]\n', 'Specific Energy', props.specEnergy/(1*hour));
            fprintf('%20s: %g [g]\n', 'packing mass (given)', ecs.packingMass/gram);


            fprintf('%-25s%20s%20s\n', '', 'theta100', 'theta0');
            fprintf('%-25s%20.5f%20.5f \n', ne, props.(ne).theta100, props.(ne).theta0);
            fprintf('%-25s%20.5f%20.5f \n', pe, props.(pe).theta100, props.(pe).theta0);
            
            % fprintf('%20s: %g [-]\n', 'N/P ratio', props.NPratio);
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
