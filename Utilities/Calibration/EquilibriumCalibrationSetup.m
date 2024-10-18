classdef EquilibriumCalibrationSetup
%% The goal of this class is to solve the following calibration problem:
% Given a discharge curve, find the set of parameters that give the better match.
% By default (calibrationCase = 1), the calibration parameters are the guest stoichiometry at discharge start and the volume fraction, for both electrodes.
%
%  Usage description :
%  1. Instantiate object using model and expdata
%     expdata is a structure with fields
%     - I    : Current used (scalar)
%     - U    : vector with voltage values
%     - time : vector with time values (same dimension as U)
%
%  2. Run one of the optimiser, that is either
%     [Xopt, hist] = runUnitBoxBFGS(ecs, X0);
%     or
%     [Xopt, info] = runIpOpt(ecs, ipopt_options);
%     (see documentation for those)
%
%  3. The optimal parameter values are given by the vector Xopt.
%     They can be converted to a more readable structure using
%     vals = assignFromX(ecs, X)
%     (see documentation for this method)
%
% See installation instruction for ipopt optimiser in README.org in this directory
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

        calibrationCase = 1
        % different calibration case depending on the parameters that are chosen. See method printVariableChoice below
        % At the moment, the following is implemented
        % case 1 (default) : The calibration parameters are theta100 and volume fraction for both electrodes
        %                    Theta0 for the positive electrode is computed from the end point of the discharge curve
        %                    Theta0 for the negative electrode is computed to match a given NP ration (default value 1.1)
        %                    When using ipopt, we add a constraint that enforces that the theta value at the end (t = totalTime) is between 0 and 1.
        % case 2           : The calibration parameters are theta100 for the negative electrode, the volume fractions for both electrodes
        % case 3           : The calibration parameters are theta100, theta0 and volume fraction for both electrodes and we add a constraint on the np-ratio (thus we use IpOpt solver)


        %% Helper structures, assigned during setup

        bounds % Vector of variable bounds on the variables, with field
               % - lower
               % - upper

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

                data_default = struct('ne_theta_max', 1, ...
                                      'np_ratio', 1.1);

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

            opt = struct('verbose', true, ...
                         'lower', 1e-3, ...
                         'upper', 1);
            opt = merge_options(opt, varargin{:});

            N = numel(ecs.getDefaultValue());

            lower = opt.lower * ones(N, 1);
            upper = opt.upper * ones(N, 1);

            ecs.bounds = struct('lower', lower, ...
                                'upper', upper);

            if opt.verbose
                fprintf('The variable bounds are\n');
                fprintf('lower\n')
                display(ecs.bounds.lower);
                fprintf('upper\n')
                display(ecs.bounds.upper);
            end

        end


        function X = getDefaultValue(ecs)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};

            switch ecs.calibrationCase

              case 1

                for ielde = 1 : numel(eldes)
                    elde     = eldes{ielde};
                    compInds = ecs.model.(elde).(co).compInds;

                    vals.(elde).theta100       = ecs.model.(elde).(co).(am).(itf).guestStoichiometry100;
                    vals.(elde).volumeFraction = ecs.model.(elde).(co).volumeFraction*ecs.model.(elde).(co).volumeFractions(compInds.(am));
                end

              case 2

                %% recall : ordering of parameters
                %
                % X(1) : theta100 anode
                % X(2) : volume fraction anode
                % X(3) : volume fraction cathode


                compInds = ecs.model.(ne).(co).compInds;
                vals.(ne).theta100       = ecs.model.(ne).(co).(am).(itf).guestStoichiometry100;
                vals.(ne).volumeFraction = ecs.model.(ne).(co).volumeFraction*ecs.model.(ne).(co).volumeFractions(compInds.(am));

                compInds = ecs.model.(pe).(co).compInds;
                vals.(pe).volumeFraction = ecs.model.(pe).(co).volumeFraction*ecs.model.(pe).(co).volumeFractions(compInds.(am));

              case 3

                % X(1) : theta100 anode
                % X(2) : theta0 anode
                % X(3) : volume fraction anode
                % X(4) : theta100 cathode
                % X(5) : theta0 cathode
                % X(6) : volume fraction cathode

                for ielde = 1 : numel(eldes)
                    elde     = eldes{ielde};
                    compInds = ecs.model.(elde).(co).compInds;

                    vals.(elde).theta100       = ecs.model.(elde).(co).(am).(itf).guestStoichiometry100;
                    vals.(elde).theta0         = ecs.model.(elde).(co).(am).(itf).guestStoichiometry0;
                    vals.(elde).volumeFraction = ecs.model.(elde).(co).volumeFraction*ecs.model.(elde).(co).volumeFractions(compInds.(am));
                end

              otherwise
                error('calibrationCase not recognized')
            end

            X = ecs.assignToX(vals);

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


        function X = assignToX(ecs, vals)

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

                X = nan(4, 1);

                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    X(2*ielde - 1) = vals.(elde).theta100;
                    X(2*ielde)     = vals.(elde).volumeFraction;
                end

              case 2

                %% recall : ordering of parameters
                %
                % X(1) : theta100 anode
                % X(2) : volume fraction anode
                % X(3) : volume fraction cathode


                X(1) = vals.(ne).theta100;
                X(2) = vals.(ne).volumeFraction;
                X(3) = vals.(pe).volumeFraction;

              case 3

                % X(1) : theta100 anode
                % X(2) : theta0 anode
                % X(3) : volume fraction anode
                % X(4) : theta100 cathode
                % X(5) : theta0 cathode
                % X(6) : volume fraction cathode

                X = nan(6, 1);

                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    X(3*ielde - 2) = vals.(elde).theta100;
                    X(3*ielde - 1) = vals.(elde).theta0;
                    X(3*ielde)     = vals.(elde).volumeFraction;
                end

              otherwise
                error('calibrationCase not recognized')
            end
        end


        function vals = setupAlphas(ecs, vals)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                vol  = sum(ecs.model.(elde).(co).G.getVolumes());
                cmax = ecs.model.(elde).(co).(am).(itf).saturationConcentration;

                vals.(elde).alpha = vals.(elde).volumeFraction*vol*cmax;
                vals.(elde).cmax  = cmax;

            end

        end


        function vals = assignFromX(ecs, X)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

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

                    vals.(elde).theta100       = X(2*ielde - 1);
                    vals.(elde).volumeFraction = X(2*ielde);

                end

              case 2

                %% recall : ordering of parameters
                %
                % X(1) : theta100 anode
                % X(2) : volume fraction anode
                % X(3) : volume fraction cathode

                vals.(ne).theta100       = X(1);
                vals.(ne).volumeFraction = X(2);
                vals.(pe).volumeFraction = X(3);

              case 3

                % X(1) : theta100 anode
                % X(2) : theta0 anode
                % X(3) : volume fraction anode
                % X(4) : theta100 cathode
                % X(5) : theta0 cathode
                % X(6) : volume fraction cathode


                for ielde = 1 : numel(eldes)

                    elde = eldes{ielde};

                    vals.(elde).theta100       = X(3*ielde - 2);
                    vals.(elde).theta0         = X(3*ielde - 1);
                    vals.(elde).volumeFraction = X(3*ielde);

                end


              otherwise

                error('calibration case not recognized');

            end

        end


        function [fexp, fcomp] = setupfunction(ecs)
        % Setup the function to compute the discharge curves, either for the experimental data or for the given set of parameter, using the model.
            fcomp = @(t, X) ecs.computeF(t, X);
            fexp  = @(t) ecs.experimentalF(t);

        end


        function theta = computeTheta(ecs, t, elde, tf, theta_tf, alpha)
        % theta_tf is lithiation at time tf*totalTime
        % returns lithiation at time t
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';

            I = ecs.expI;
            T = ecs.totalTime;

            switch elde
              case ne
                sgn = -1;
              case pe
                sgn = 1;
            end

            theta = theta_tf + (t - tf*T)*((sgn*I) ./ (ecs.F*alpha));

        end


        function [f, fpe, fne] = computeF(ecs, t, X)

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

            T     = ecs.Temperature;

            vals = ecs.updateThetas(X, 'includeTheta0', false);

            theta100 = vals.(pe).theta100;
            alpha    = vals.(pe).alpha;

            theta = ecs.computeTheta(t, pe, 0, theta100, alpha);
            cmax = vals.(pe).cmax;
            fpe = ecs.model.(pe).(co).(am).(itf).computeOCPFunc(theta*cmax, T, cmax);

            theta100 = vals.(ne).theta100;
            alpha    = vals.(ne).alpha;

            theta = ecs.computeTheta(t, ne, 0, theta100, alpha);
            cmax = vals.(ne).cmax;
            fne = ecs.model.(ne).(co).(am).(itf).computeOCPFunc(theta*cmax, T, cmax);

            f = fpe - fne;

        end


        function f = experimentalF(ecs, t)

            time = ecs.exptime; % in second
            U    = ecs.expU;

            f = interp1(time, U, t, 'linear', 'extrap');

        end


        function [z, dz] = objective(ecs, X)

            t  = ecs.exptime;

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

            vals = ecs.updateThetas(X, 'includeTheta0', true);

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            np_ratio = (vals.(ne).alpha*(vals.(ne).theta0 - vals.(ne).theta100))./(vals.(pe).alpha*(vals.(pe).theta100 - vals.(pe).theta0));

        end


        function y = constraints_func(ecs, X)
        % for ipopt api

            switch ecs.calibrationCase

              case 1

                % We enforce that theta at total time is between 0 and 1

                I = ecs.expI;
                T = ecs.totalTime;

                vals = ecs.assignFromX(X);
                vals = ecs.setupAlphas(vals);

                ne  = 'NegativeElectrode';
                pe  = 'PositiveElectrode';

                eldes = {ne, pe};

                y = {};

                for ielde = 1 : numel(eldes)

                    elde = eldes{ielde};

                    alpha    = vals.(elde).alpha;
                    theta100 = vals.(elde).theta100;

                    switch elde
                      case ne
                        sgn = -1;
                      case pe
                        sgn = 1;
                    end

                    theta = theta100 + T*((sgn*I)./(ecs.F*alpha));

                    y{end + 1} = theta;
                    y{end + 1} = 1 - theta;

                end

                y = vertcat(y{:});

              case 2

                error('not implemented');

              case 3

                data = ecs.calibrationParameters;
                np_ratio = ecs.computeNPratio(X);

                y = np_ratio./data.np_ratio - 1;

              otherwise

                error('ecs.calibrationCase not recognized');

            end
        end


        function dy = jacobian_constraints_func(ecs, X)
        % for ipopt api

            X = initVariablesADI(X);
            y = ecs.constraints_func(X);
            dy = y.jac{1};

        end


        function y = jacobian_constraints_structure_func(ecs)
        % we do not bother here and setup as full

            switch ecs.calibrationCase

              case 1

                y = sparse(ones(4, 4));

              case 2

                error('not implemented');

              case 3

                y = sparse(ones(1, 6));

              otherwise

                error('ecs.calibrationCase not recognized');

            end
        end


        function vals = updateThetas(ecs, X, varargin)

            opt = struct('includeTheta0', false);
            opt = merge_options(opt, varargin{:});

            totalTime = ecs.totalTime;

            vals = ecs.assignFromX(X);
            vals = ecs.setupAlphas(vals);

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            eldes = {ne, pe};

            switch ecs.calibrationCase

              case 1

                if opt.includeTheta0

                    % We compute theta0 in cathode as the lithiation at end of discharge in
                    vals.(pe).theta0 = ecs.computeTheta(totalTime, pe, 0, vals.(pe).theta100, vals.(pe).alpha);

                    % vals.(ne).theta0 = ecs.computeTheta(totalTime, ne, 0, vals.(ne).theta100, vals.(ne).alpha);

                    data = ecs.calibrationParameters;
                    np_ratio = data.np_ratio;

                    vals.(ne).theta0 = vals.(ne).theta100 - np_ratio*(vals.(pe).alpha/vals.(ne).alpha)*(vals.(pe).theta0 - vals.(pe).theta100);

                else

                    % nothing to do

                end

              case 2

                data = ecs.calibrationParameters;
                vals.(pe).theta100 = data.pe_theta;

                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    vals.(elde).theta0 = ecs.computeTheta(totalTime, elde, 0, vals.(elde).theta100, vals.(elde).alpha);
                end

              case 3

                % nothing to do

              otherwise

                error('calibration case not recognized');

            end

        end


        function vals = computeCapacitiesFromVals(ecs, vals)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                cmax  = ecs.model.(elde).(co).(am).(itf).saturationConcentration;

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

                cap = (cM - cm)/cmax*alpha*ecs.F;

                vals.(elde).cmax = cmax;
                vals.(elde).cap  = cap;

            end

        end


        function vals = computeCapacities(ecs, X)

            vals = ecs.updateThetas(X, 'includeTheta0', true);
            vals = ecs.computeCapacitiesFromVals(vals);

        end


        function props = computeProperties(ecs, X)

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

                props.(elde).dischargeFunc = @(s) ecs.model.(elde).(co).(am).(itf).computeOCPFunc(c, T, cmax);

                energies.(elde) = cap*smax/N*sum(props.(elde).dischargeFunc(s));

                props.(elde).U = ecs.model.(elde).(co).(am).(itf).computeOCPFunc(c0, T, cmax);

            end

            props.energy = energies.(pe) - energies.(ne);

            mass = ecs.getMass();
            props.specEnergy = props.energy/mass;

            % props.NPratio = props.(ne).cap/props.(pe).cap;

            props.U = props.(pe).U - props.(ne).U;

        end


        function mass = getMass(ecs)

            if isempty(ecs.packingMass)
                fprintf('Packing Mass is not given, we set it equal to zero\n');
                ecs.packingMass = 0;
            end

            mass = computeCellMass(ecs.model, 'packingMass', ecs.packingMass);

        end


        function [Xopt, hist] = runUnitBoxBFGS(ecs, varargin)

            opt = struct('X0', ecs.getDefaultValue());
            [opt, extra] = merge_options(opt, varargin{:});

            f = @(X) ecs.objective(X);

            n = size(ecs.bounds.lower, 1);
            linIneq.A = [-eye(n); eye(n)];
            linIneq.b = [-ecs.bounds.lower; ecs.bounds.upper];

            params = {'objChangeTol'    , 1e-12, ...
                      'maximize'        , false, ...
                      'maxit'           , 1000 , ...
                      'maxInitialUpdate', 1e-6 , ...
                      'enforceFeasible' , true , ...
                      'lineSearchMaxIt' , 10   , ...
                      'wolfe2'          , 0.99 , ...
                      'linIneq'         , linIneq};

            % NB: will prefer options in extra over params
            [~, Xopt, hist] = unitBoxBFGS(opt.X0, f, params{:}, extra{:});

        end


        function [Xopt, info] = runIpOpt(ecs, ipopt_options)

            X0 = ecs.getDefaultValue();

            funcs.objective         = @(x) ecs.objective_func(x);
            funcs.gradient          = @(x) ecs.gradient_func(x);
            funcs.constraints       = @(x) ecs.constraints_func(x);
            funcs.jacobian          = @(x) ecs.jacobian_constraints_func(x);
            funcs.jacobianstructure = @()  ecs.jacobian_constraints_structure_func();

            options.lb = ecs.bounds.lower;
            options.ub = ecs.bounds.upper;

            switch ecs.calibrationCase

              case 1

                options.cl = 0*ones(4, 1);
                options.cu = inf*ones(4, 1);

              case 2

                error('not implemented');

              case 3

                options.cl = 0;
                options.cu = 0;

              otherwise

                error('ecs.calibrationCase not recognized');

            end

            if nargin > 1
                options.ipopt = ipopt_options;
            end
            options.ipopt.hessian_approximation = 'limited-memory';

            [Xopt, info] = ipopt(X0, funcs, options);

        end


        function printParameters(ecs, X)

            vals = ecs.assignFromX(X);

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';

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


        function json = export(ecs, X)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            itf = 'Interface';
            am  = 'ActiveMaterial';
            co  = 'Coating';

            props = ecs.computeProperties(X);

            eldes = {ne, pe};

            for ielde = 1:numel(eldes)

                elde = eldes{ielde};

                json.(elde).(co).(am).(itf).guestStoichiometry0   = props.(elde).theta0;
                json.(elde).(co).(am).(itf).guestStoichiometry100 = props.(elde).theta100;
                json.(elde).(co).volumeFraction                   = props.(elde).volumeFraction;

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
