classdef EquilibriumCalibrationSetup
%% The goal of this class is to solve the following calibration problem:
% Given a discharge curve, find the set of parameters that give the better match.
% By default (calibrationCase = 1), the calibration parameters are the guest stoichiometry at discharge start and the total amount Lithium, for both electrodes.
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
        % case 1 (default) : The calibration parameters are guestStoichiometry100 and total amount Lithium for both electrodes
        %                    guestStoichiometry0 for the positive electrode is computed from the end point of the discharge curve
        %                    guestStoichiometry0 for the negative electrode is computed to match a given NP ration (default value 1.1)
        %                    When using ipopt, we add a constraint that enforces that the theta value at the end (t = totalTime) is between 0 and 1.
        % case 2           : The calibration parameters are guestStoichiometry100 for the negative electrode, the total amount Lithium for both electrodes
        % case 3           : The calibration parameters are guestStoichiometry100, guestStoichiometry0 and total amount Lithium for both electrodes and we add a constraint on the np-ratio (thus we use IpOpt solver)

        lowerCutoffVoltage % value of the lower cutoff voltage that is used to compute guestStoichiometry0. if not given the value is
                           % computed from expdata

        totalAmountVariableChoice  % variable that is chosen to adjust the total amount of lithium in the electrode. The choices are
                                   % - 'volumeFraction'          : the volume fraction of the active material
                                   % - 'saturationConcentration' : the saturation concentration of the electrode

        %% Helper structures, assigned during setup

        bounds % Vector of variable bounds on the variables, with field
               % - lower
               % - upper

        calibrationParameters

        vals0 % initial values of the calibration parameters
        X0    % initial values of the calibration parameters in the form of a vector

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

            if isempty(ecs.lowerCutoffVoltage)
                ecs.lowerCutoffVoltage = min(ecs.expU);
            end

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

            ecs.vals0 = ecs.getDefaultValue();
            ecs.X0    = ecs.assignToX(ecs.vals0);

            ecs = ecs.setupDefaultVariableBounds('verbose', opt.verbose);

        end


        function ecs = setupDefaultVariableBounds(ecs, varargin)

            opt = struct('verbose'      , true, ...
                         'relativeLower', 1e-1, ...
                         'relativeUpper', 1e1);

            opt = merge_options(opt , varargin{:});

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};

            vals0 = ecs.vals0;


            % lower bound

            vals = vals0;
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                vals.(elde).guestStoichiometry100 = 0;
                vals.(elde).guestStoichiometry0   = 0;
                vals.(elde).totalAmount         = vals0.(elde).totalAmount*opt.relativeLower;

            end

            lower = ecs.assignToX(vals);

            % upper bound

            vals = vals0;
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                vals.(elde).guestStoichiometry100 = 1;
                vals.(elde).guestStoichiometry0   = 1;
                vals.(elde).totalAmount         = vals0.(elde).totalAmount*opt.relativeUpper;

            end

            upper = ecs.assignToX(vals);

            % assign to ecs

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


        function vals0 = getDefaultValue(ecs)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};

            switch ecs.calibrationCase

              case 1

                for ielde = 1 : numel(eldes)

                    elde    = eldes{ielde};
                    coating = ecs.model.(elde).(co);
                    amInd   = coating.compInds.(am);

                    vol = sum(coating.G.getVolumes());

                    vals0.(elde).guestStoichiometry100   = coating.(am).(itf).guestStoichiometry100;
                    vals0.(elde).guestStoichiometry0     = coating.(am).(itf).guestStoichiometry0; % just for reference, not used in the calibration
                    vals0.(elde).saturationConcentration = coating.(am).(itf).saturationConcentration; % used in computation of OCP.
                    vals0.(elde).totalAmount             = vol*coating.volumeFraction*coating.volumeFractions(amInd)*coating.(am).(itf).saturationConcentration;
                end

                vals0.np_ratio = (vals0.(ne).totalAmount*(vals0.(ne).guestStoichiometry0 - vals0.(ne).guestStoichiometry100))./(vals0.(pe).totalAmount*(vals0.(pe).guestStoichiometry100 - vals0.(pe).guestStoichiometry0));

              case 2

                %% recall : ordering of parameters
                %
                % X(1) : guestStoichiometry100 anode
                % X(2) : total amount Lithium anode
                % X(3) : total amount Lithium cathode


                coating = ecs.model.(ne).(co);
                amInd = coating.compInds.(am);

                vol = sum(coating.G.getVolumes());
                vals0.(ne).guestStoichiometry100   = coating.(am).(itf).guestStoichiometry100;
                vals0.(ne).guestStoichiometry0     = coating.(am).(itf).guestStoichiometry0; % just for reference, not used in the calibration
                vals0.(ne).saturationConcentration = coating.(am).(itf).saturationConcentration; % used in computation of OCP.
                vals0.(ne).totalAmount             = vol.*coating.volumeFraction*coating.volumeFractions(amInd)*coating.(am).(itf).saturationConcentration;

                coating = ecs.model.(pe).(co);
                amInd = coating.compInds.(am);

                vol = sum(coating.G.getVolumes());
                vals0.(pe).saturationConcentration = coating.(am).(itf).saturationConcentration; % used in computation of OCP.
                vals0.(pe).totalAmount = vol.*coating.volumeFraction*coating.volumeFractions(amInd)*coating.(am).(itf).saturationConcentration;

              case 3

                % X(1) : guestStoichiometry100 anode
                % X(2) : guestStoichiometry0 anode
                % X(3) : total amount Lithium anode
                % X(4) : guestStoichiometry100 cathode
                % X(5) : guestStoichiometry0 cathode
                % X(6) : total amount Lithium cathode

                for ielde = 1 : numel(eldes)

                    elde    = eldes{ielde};
                    coating = ecs.model.(elde).(co);
                    amInd   = coating.compInds.(am);

                    vol = sum(coating.G.getVolumes());
                    vals0.(elde).guestStoichiometry100   = coating.(am).(itf).guestStoichiometry100;
                    vals0.(elde).guestStoichiometry0     = coating.(am).(itf).guestStoichiometry0;
                    vals0.(elde).saturationConcentration = coating.(am).(itf).saturationConcentration; % used in computation of OCP.
                    vals0.(elde).totalAmount             = vol*coating.volumeFraction*coating.volumeFractions(amInd)*coating.(am).(itf).saturationConcentration;

                end

              otherwise

                error('calibrationCase not recognized')

            end

        end


        function printVariableChoice(ecs)

            switch ecs.calibrationCase
              case 1
                fprintf('\nThe calibration parameters are guestStoichiometry100 and total amount Lithium for both electrodes\n');
              case 2
                fprintf('\nThe calibration parameters are guestStoichiometry100 for the pegative electrode\n');
                fprintf('The total amount Lithium for both electrodes\n');
              case 3
                fprintf('\nThe calibration parameters are guestStoichiometry100, guestStoichiometry0, and total amount Lithium for both electrodes\n');
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
                % X(1) : guestStoichiometry100 anode
                % X(2) : total amount Lithium anode
                % X(3) : guestStoichiometry100 cathode
                % X(4) : total amount Lithium cathode

                X = nan(4, 1);

                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    X(2*ielde - 1) = vals.(elde).guestStoichiometry100;
                    X(2*ielde)     = vals.(elde).totalAmount;
                end

              case 2

                %% recall : ordering of parameters
                %
                % X(1) : guestStoichiometry100 anode
                % X(2) : total amount Lithium anode
                % X(3) : total amount Lithium cathode


                X(1) = vals.(ne).guestStoichiometry100;
                X(2) = vals.(ne).totalAmount;
                X(3) = vals.(pe).totalAmount;

              case 3

                % X(1) : guestStoichiometry100 anode
                % X(2) : guestStoichiometry0 anode
                % X(3) : total amount Lithium anode
                % X(4) : guestStoichiometry100 cathode
                % X(5) : guestStoichiometry0 cathode
                % X(6) : total amount Lithium cathode

                X = nan(6, 1);

                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    X(3*ielde - 2) = vals.(elde).guestStoichiometry100;
                    X(3*ielde - 1) = vals.(elde).guestStoichiometry0;
                    X(3*ielde)     = vals.(elde).totalAmount;
                end

              otherwise
                error('calibrationCase not recognized')
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
                % X(1) : guestStoichiometry100 anode
                % X(2) : total amount Lithium anode
                % X(3) : guestStoichiometry100 cathode
                % X(4) : total amount Lithium cathode

                for ielde = 1 : numel(eldes)

                    elde = eldes{ielde};

                    vals.(elde).guestStoichiometry100 = X(2*ielde - 1);
                    vals.(elde).totalAmount           = X(2*ielde);

                end

              case 2

                %% recall : ordering of parameters
                %
                % X(1) : guestStoichiometry100 anode
                % X(2) : total amount Lithium anode
                % X(3) : total amount Lithium cathode

                vals.(ne).guestStoichiometry100 = X(1);
                vals.(ne).totalAmount           = X(2);
                vals.(pe).totalAmount           = X(3);

              case 3

                % X(1) : guestStoichiometry100 anode
                % X(2) : guestStoichiometry0 anode
                % X(3) : total amount Lithium anode
                % X(4) : guestStoichiometry100 cathode
                % X(5) : guestStoichiometry0 cathode
                % X(6) : total amount Lithium cathode


                for ielde = 1 : numel(eldes)

                    elde = eldes{ielde};

                    vals.(elde).guestStoichiometry100 = X(3*ielde - 2);
                    vals.(elde).guestStoichiometry0   = X(3*ielde - 1);
                    vals.(elde).totalAmount           = X(3*ielde);

                end


              otherwise

                error('calibration case not recognized');

            end

        end


        function [fexp, fcomp] = setupfunction(ecs)
        % Setup the function to compute the discharge curves, either for the experimental data or for the given set of
        % parameter, using the model.

            fcomp = @(t, X) ecs.computeF(t, X);
            fexp  = @(t) ecs.experimentalF(t);

        end


        function theta = computeTheta(ecs, t, elde, tf, theta_tf, totalAmount)
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

            theta = theta_tf + (t - tf*T)*((sgn*I) ./ (ecs.F*totalAmount));

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

            T = ecs.Temperature;

            vals  = ecs.assignFromX(X);
            vals0 = ecs.vals0;

            guestStoichiometry100 = vals.(pe).guestStoichiometry100;
            totalAmount           = vals.(pe).totalAmount;

            theta = ecs.computeTheta(t, pe, 0, guestStoichiometry100, totalAmount);
            cmax  = vals0.(pe).saturationConcentration;
            fpe = ecs.model.(pe).(co).(am).(itf).computeOCPFunc(theta*cmax, T, cmax);

            guestStoichiometry100 = vals.(ne).guestStoichiometry100;
            totalAmount           = vals.(ne).totalAmount;

            theta = ecs.computeTheta(t, ne, 0, guestStoichiometry100, totalAmount);
            cmax  = vals0.(ne).saturationConcentration;
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

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            vals = ecs.updateGuestStoichiometries(X, 'includeGuestStoichiometry0', true);

            np_ratio = (vals.(ne).totalAmount*(vals.(ne).guestStoichiometry0 - vals.(ne).guestStoichiometry100))./(vals.(pe).totalAmount*(vals.(pe).guestStoichiometry100 - vals.(pe).guestStoichiometry0));

        end


        function y = constraints_func(ecs, X)
        % for ipopt api

            switch ecs.calibrationCase

              case 1

                % We enforce that theta at total time is between 0 and 1

                I = ecs.expI;
                T = ecs.totalTime;

                vals = ecs.assignFromX(X);

                ne  = 'NegativeElectrode';
                pe  = 'PositiveElectrode';

                eldes = {ne, pe};

                y = {};

                for ielde = 1 : numel(eldes)

                    elde = eldes{ielde};

                    totalAmount           = vals.(elde).totalAmount;
                    guestStoichiometry100 = vals.(elde).guestStoichiometry100;

                    switch elde
                      case ne
                        sgn = -1;
                      case pe
                        sgn = 1;
                    end

                    theta = guestStoichiometry100 + T*((sgn*I)./(ecs.F*totalAmount));

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


        function vals = updateGuestStoichiometries(ecs, X, varargin)

            opt = struct('includeGuestStoichiometry0', false);
            opt = merge_options(opt, varargin{:});

            totalTime = ecs.totalTime;

            vals = ecs.assignFromX(X);

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            eldes = {ne, pe};

            switch ecs.calibrationCase

              case 1

                if opt.includeGuestStoichiometry0

                    [~, fcomp] = ecs.setupfunction();
                    tend = totalTime;
                    done = false;
                    while ~done
                        Uend = fcomp(tend, X);
                        if Uend < ecs.lowerCutoffVoltage
                            done = true;
                        else
                            tend = 1.5*tend;
                        end
                    end
                    t = linspace(0, tend, 10000);

                    tend = t(fcomp(t, X) > ecs.lowerCutoffVoltage);
                    tend = tend(end);

                    vals.(pe).guestStoichiometry0 = ecs.computeTheta(tend, pe, 0, vals.(pe).guestStoichiometry100, vals.(pe).totalAmount);
                    if vals.(pe).guestStoichiometry0 > 1
                        warning('Set guestStoichiometry0 of positive electrode to 1 (from %g)', vals.(pe).guestStoichiometry0);
                        vals.(pe).guestStoichiometry0 = 1;
                    end

                    data = ecs.calibrationParameters;
                    np_ratio = data.np_ratio;

                    vals.(ne).guestStoichiometry0 = vals.(ne).guestStoichiometry100 - np_ratio*(vals.(pe).totalAmount/vals.(ne).totalAmount)*(vals.(pe).guestStoichiometry0 - vals.(pe).guestStoichiometry100);
                    if vals.(ne).guestStoichiometry0 < 0
                        warning('Set guestStoichiometry0 of negative electrode to 0 (from %g)', vals.(ne).guestStoichiometry0);
                        vals.(ne).guestStoichiometry0 = 0;
                        vals.np_ratio = vals.(ne).guestStoichiometry100/((vals.(pe).totalAmount/vals.(ne).totalAmount)*(vals.(pe).guestStoichiometry0 - vals.(pe).guestStoichiometry100));
                    end

                else

                    % nothing to do

                end

              case 2

                data = ecs.calibrationParameters;
                vals.(pe).guestStoichiometry100 = data.pe_theta;

                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    vals.(elde).guestStoichiometry0 = ecs.computeTheta(totalTime, elde, 0, vals.(elde).guestStoichiometry100, vals.(elde).totalAmount);
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

                guestStoichiometry0   = vals.(elde).guestStoichiometry0;
                guestStoichiometry100 = vals.(elde).guestStoichiometry100;
                totalAmount    = vals.(elde).totalAmount;

                switch elde
                  case ne
                    cM = guestStoichiometry100*cmax;
                    cm = guestStoichiometry0*cmax;
                  case pe
                    cM = guestStoichiometry0*cmax;
                    cm = guestStoichiometry100*cmax;
                end

                cap = (cM - cm)/cmax*totalAmount*ecs.F;

                vals.(elde).saturationConcentration = cmax;
                vals.(elde).cap  = cap;

            end

        end

        function vals = computeCapacities(ecs, X)

            vals = ecs.updateGuestStoichiometries(X, 'includeGuestStoichiometry0', true);
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
                cmax = props.(elde).saturationConcentration;
                c0   = props.(elde).guestStoichiometry100*cmax;
                cT   = props.(elde).guestStoichiometry0*cmax;
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

            opt = struct('X0', ecs.X0);
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

            X0 = ecs.X0;

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
                fprintf('%-25s%25s%25s\n', '', 'guestStoichiometry100', 'total amount Lithium');
                thetastr = sprintf('%6.5f', vals.(ne).guestStoichiometry100);
                fprintf('%-25s%25s%25.5f \n', ne, thetastr, vals.(ne).totalAmount);
                thetastr = sprintf('%6.5f', vals.(pe).guestStoichiometry100);
                fprintf('%-25s%25s%25.5f \n', pe, thetastr, vals.(pe).totalAmount);
              case 2
                fprintf('%-25s%20s%20s\n', '', 'guestStoichiometry100', 'total amount Lithium');
                thetastr = sprintf('%6.5f', vals.(ne).guestStoichiometry100);
                fprintf('%-25s%20s%20.5f \n', ne, thetastr, vals.(ne).totalAmount);
                fprintf('%-25s%20s%20.5f \n', pe, '', vals.(pe).totalAmount);
              case 3
                fprintf('%-25s%20s%20s%20s\n', '', 'guestStoichiometry100', 'guestStoichiometry0', 'total amount Lithium');
                fprintf('%-25s%20.5f%20.5f%20.5f \n', ne, vals.(ne).guestStoichiometry100, vals.(ne).guestStoichiometry0, vals.(ne).totalAmount);
                fprintf('%-25s%20.5f%20.5f%20.5f \n', pe, vals.(pe).guestStoichiometry100, vals.(pe).guestStoichiometry0, vals.(pe).totalAmount);
              otherwise
                error('calibration case not recognized');
            end

        end

        function printProperties(ecs, X)

            props = ecs.computeProperties(X);

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';

            fprintf('%20s: %g [V]\n'    , 'initial voltage'     , props.U);
            fprintf('%20s: %g [Ah]\n'   , 'Capacity'            , props.cap/(1*hour));
            fprintf('%20s: %g [Wh/kg]\n', 'Specific Energy'     , props.specEnergy/(1*hour));
            fprintf('%20s: %g [g]\n'    , 'packing mass (given)', ecs.packingMass/gram);


            fprintf('%-25s%20s%20s\n'     , '', 'guestStoichiometry100'         , 'guestStoichiometry0');
            fprintf('%-25s%20.5f%20.5f \n', ne, props.(ne).guestStoichiometry100, props.(ne).guestStoichiometry0);
            fprintf('%-25s%20.5f%20.5f \n', pe, props.(pe).guestStoichiometry100, props.(pe).guestStoichiometry0);

            % fprintf('%20s: %g [-]\n', 'N/P ratio', props.NPratio);
        end

        function jsonstruct = exportParameters(ecs, X, varargin)

            opt = struct('guestStoichiometries0', []);
            opt = merge_options(opt, varargin{:});

            assert(~isempty(ecs.totalAmountVariableChoice), 'The total amount of Lithium is given by the product of the volume, the volume fraction, the active material volume fraction and the saturation concentration. To setup the model, we have to choose one of these variables. The choice is made by setting the field totalAmountVariableChoice in the EquilibriumCalibrationSetup structure. You can either set it to  ''volumeFraction'' or ''saturationConcentration''');


            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};

            vals = ecs.updateGuestStoichiometries(X, 'includeGuestStoichiometry0', true);

            vals0 = ecs.vals0;

            jsonstruct = struct();

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                guestStoichiometry100 = vals.(elde).guestStoichiometry100;

                if isAssigned(opt.guestStoichiometries0, {elde})
                    guestStoichiometry0 = opt.guestStoichiometries0.(elde);
                else
                    guestStoichiometry0 = vals.(elde).guestStoichiometry0;
                end

                totalAmount  = vals.(elde).totalAmount;
                totalAmount0 = vals0.(elde).totalAmount;

                jsonstruct = setJsonStructField(jsonstruct, {elde, co, am, itf, 'guestStoichiometry100'}, guestStoichiometry100);
                jsonstruct = setJsonStructField(jsonstruct, {elde, co, am, itf, 'guestStoichiometry0'}, guestStoichiometry0);

                switch ecs.totalAmountVariableChoice

                  case 'volumeFraction'

                    jsonstruct = setJsonStructField(jsonstruct, {elde, co, 'volumeFraction'}, ...
                                                    (totalAmount/totalAmount0) * ecs.model.(elde).(co).volumeFraction);

                  case 'saturationConcentration'

                    jsonstruct = setJsonStructField(jsonstruct, {elde, co, am, itf, 'saturationConcentration'}, ...
                                                    (totalAmount/totalAmount0)*ecs.model.(elde).(co).(am).(itf).saturationConcentration);

                  otherwise

                    error('Unrecognized totalAmountVariableChoice %s', ecs.totalAmountVariableChoice);

                end

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
