classdef EquilibriumCalibrationSetup

    properties
        
        model % Battery model

        F % Faraday Constant
        Temperature = 298.15% Temperature (used in OCP curves)

        % Discharge curve measured experementally at near-equilibrium condition
        exptime % Time [s]
        expU    % Voltage [V]
        expI    % Current [A]
        totalTime % total time (set as exptime(end))
        
        packingMass = 61*gram % mass of packing 

        calibrationCase = 1

        calibrationCaseParameters 
        
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

            % setup some default
            data = struct('ne_theta_max', 1);
            calibrationCaseParameters{1} = data;
            data = struct('pe_theta', 0.99, ...
                          'ne_theta_max', 1);
            calibrationCaseParameters{2} = data;

            ecs.calibrationCaseParameters = calibrationCaseParameters;

        end

        function data = getCaseParameters(ecs)

            data = ecs.calibrationCaseParameters{ecs.calibrationCase};
            
        end

        function ecs = setCaseParameters(ecs, data)

            ecs.calibrationCaseParameters{ecs.calibrationCase} = data;
            
        end

        function linIneq = setupIneqConstraints(ecs)

            data = ecs.getCaseParameters();
            
            switch ecs.calibrationCase
                
              case 1

                A = -eye(4);
                b = -0.1*[1, 1, 1, 1]';

                A = [A; 1, 0, 0, 0];
                b = [b; data.ne_theta_max];

              case 2
                
                A = -eye(3);
                b = -0.1*[1, 1, 1]';
                
                A = [A; 1, 0, 0];
                b = [b; data.ne_theta_max];
                
            end

            linIneq.A = A;
            linIneq.b = b;
            
        end

        function params = setupOptimParams(ecs)
            
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
        
        function ecs = setPackingMass(ecs, packingMass)

            model = ecs.model;
            
            [mass, masses] = computeCellMass(model, 'packingMass', packingMass);

            ecs.packingMass = packingMass;
            ecs.mass        = mass;
            ecs.masses      = masses;
            
        end

        function X = getDefaultValue(ecs)
        %% recall : ordering of parameters
        %
        % X(1) : theta100 anode
        % X(2) : alpha anode (alpha = V*volumeFraction*cmax)
        % X(3) : theta100 cathode
        % X(4) : alpha cathode (alpha = V*volumeFraction*cmax)

            
        % We define some shorthand names for simplicity.
            model = ecs.model;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            switch ecs.calibrationCase

              case 1
                
                X = nan(4, 1);

                X(1) = model.(ne).(am).(itf).theta100;
                vf   = unique(model.(ne).(am).volumeFraction)*model.(ne).(am).activeMaterialFraction;
                X(2) = vf;
                
                X(3) = model.(pe).(am).(itf).theta100;
                vf   = unique(model.(pe).(am).volumeFraction)*model.(pe).(am).activeMaterialFraction;
                X(4) = vf;

              case 2
                
                X = nan(3, 1);

                X(1) = model.(ne).(am).(itf).theta100;
                vf   = unique(model.(ne).(am).volumeFraction)*model.(ne).(am).activeMaterialFraction;
                X(2) = vf;
                
                vf   = unique(model.(pe).(am).volumeFraction)*model.(pe).(am).activeMaterialFraction;
                X(3) = vf;
                
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
                fprintf('\nThe volume fractions for both electrodes\n');
              otherwise
                error('calibrationCase not recognized');
            end
            
               
        end
        
        function vals = getPhysicalValues(ecs, X)

            model = ecs.model;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};
            
            switch ecs.calibrationCase

              case 1

                for ielde = 1 : numel(eldes)

                    elde = eldes{ielde};
                    
                    vol  = sum(model.(elde).(am).G.cells.volumes);
                    cmax = model.(elde).(am).(itf).cmax;

                    vals.(elde).tf             = 0;
                    vals.(elde).theta          = X(2*ielde - 1);
                    vals.(elde).alpha          = X(2*ielde)*vol*cmax;
                    vals.(elde).volumeFraction = X(2*ielde);
                    
                end

              case 2
                
                vol  = sum(model.(ne).(am).G.cells.volumes);
                cmax = model.(ne).(am).(itf).cmax;
                
                vals.(ne).tf             = 0;
                vals.(ne).theta          = X(1);
                vals.(ne).volumeFraction = X(2);
                vals.(ne).alpha          = X(2)*vol*cmax;

                vol  = sum(model.(pe).(am).G.cells.volumes);
                cmax = model.(pe).(am).(itf).cmax;
                
                vals.(pe).volumeFraction = X(3);
                vals.(pe).alpha          = X(3)*vol*cmax;
                
                data = ecs.getCaseParameters();
                vals.(pe).tf             = 1;
                vals.(pe).theta          = data.pe_theta;
                
            end
                
        end
        
        
        function [fexp, fcomp] = setupfunction(ecs)

            fcomp = @(t, X) ecs.computeF(t, X);
            fexp  = @(t) ecs.experimentalF(t);
            
        end


        function theta = conc(ecs, t, elde, tf, theta_tf, alpha)
        % theta is lithiation at time tf*totalTimr
        % returns lithiation
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

        %% ordering of parameters
        %
        % X(1) : theta100 cathode
        % X(2) : alpha cathode (alpha = V*volumeFraction*cmax)
        % X(3) : theta100 anode
        % X(4) : alpha anode (alpha = V*volumeFraction*cmax)


            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = "SolidDiffusion";
            thermal = 'ThermalModel';
            ctrl    = 'Control';
            
            I     = ecs.expI;
            model = ecs.model;
            T     = ecs.Temperature;

            vals = ecs.getPhysicalValues(X);
            
            theta_tf = vals.(pe).theta;
            tf       = vals.(pe).tf;
            alpha    = vals.(pe).alpha;

            theta = ecs.conc(t, pe, tf, theta_tf, alpha);
            f = model.(pe).(am).(itf).computeOCPFunc(theta, T, 1);
            
            theta_tf = vals.(ne).theta;
            tf       = vals.(ne).tf;
            alpha    = vals.(ne).alpha;
            
            theta = ecs.conc(t, ne, tf, theta_tf, alpha);
            f = f - model.(ne).(am).(itf).computeOCPFunc(theta, T, 1);

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

        function props = computeProperties(ecs, X)

            model = ecs.model;
            F         = ecs.F;
            T         = ecs.Temperature;
            totalTime = ecs.totalTime;
            
            props = ecs.getPhysicalValues(X);
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};
                
                cmax  = model.(elde).(am).(itf).cmax;
                tf    = props.(elde).tf;
                theta = props.(elde).theta;
                alpha = props.(elde).alpha;
                
                theta0   = ecs.conc(totalTime, elde, tf, theta, alpha);
                theta100 = ecs.conc(0, elde, tf, theta, alpha);
                
                switch elde
                  case ne
                    cM = theta100*cmax;
                    cm = theta0*cmax;
                  case pe
                    cM = theta0*cmax;
                    cm = theta100*cmax;
                end

                cap = (cM - cm)/cmax*alpha*F;
                
                props.(elde).theta0   = theta0;
                props.(elde).theta100 = theta100;
                props.(elde).cmax     = cmax;
                props.(elde).cap      = cap;
                
            end
            
            props.cap = min(props.(ne).cap, props.(pe).cap);
            
            props.(ne).r = props.cap/props.(ne).cap;
            props.(pe).r = props.cap/props.(pe).cap;

            N = 1000; % discretization parameter
                        
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                
                smax = props.(elde).r;
                cmax = props.(elde).cmax;
                c0   = props.(elde).theta100*cmax;
                cT   = props.(elde).theta0*cmax;
                cap  = props.(elde).cap;

                s = smax.*linspace(0, 1, N + 1)';
                c = (1 - s).*c0 + s.*cT;
                f = model.(elde).(am).(itf).computeOCPFunc(c(1 : end - 1), 298, cmax);

                props.(elde).dischargeFunc = @(s) model.(elde).(am).(itf).computeOCPFunc((1 - s).*c0 + s.*cT, T, cmax);
                
                energies.(elde) = cap*smax/N*sum(props.(elde).dischargeFunc(s));

                props.(elde).U = model.(elde).(am).(itf).computeOCPFunc(c0, T, cmax);
                
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
            
            mass = computeCellMass(model, 'packingMass', packingMass);

        end
        
        function printParameters(ecs, X)

            vals = ecs.getPhysicalValues(X);
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';

            eldes = {ne, pe};


            switch ecs.calibrationCase
              case 1
                fprintf('%-25s%20s%20s\n', '', 'theta', 'volume fraction');
                thetastr = sprintf('%6.5f (SOC = %3.0f%%)', vals.(ne).theta , vals.(ne).tf*100);
                fprintf('%-25s%20s%20.5f \n', ne, thetastr, vals.(ne).volumeFraction);
                thetastr = sprintf('%6.5f (SOC = %3.0f%%)', vals.(pe).theta , vals.(pe).tf*100);
                fprintf('%-25s%20s%20.5f \n', pe, thetastr, vals.(pe).volumeFraction);
              case 2
                fprintf('%-25s%20s%20s\n', '', 'theta', 'volume fraction');
                thetastr = sprintf('%6.5f (SOC = %3.0f%%)', vals.(ne).theta , vals.(ne).tf*100);
                fprintf('%-25s%20s%20.5f \n', ne, thetastr, vals.(ne).volumeFraction);
                fprintf('%-25s%20s%20.5f \n', pe, '', vals.(pe).volumeFraction);
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
