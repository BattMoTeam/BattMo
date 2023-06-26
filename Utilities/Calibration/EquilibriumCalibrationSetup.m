classdef EquilibriumCalibrationSetup

    properties
        
        F % Faraday Constant
        exptime
        expU
        totalTime
        I
        model
        T % Temperature
        
        packingMass
        mass
        masses
        
    end

    methods

        function ecs = EquilibriumCalibrationSetup(model, expdata)

            con = PhysicalConstants();
            
            ecs.F         = con.F;
            ecs.exptime   = expdata.time;
            ecs.expU      = expdata.U;
            ecs.I         = expdata.I;
            ecs.totalTime = ecs.exptime(end);
            ecs.model     = model;

            ecs.T = 298.15; % default value
            
            packingMass = 61*gram; % default value

            ecs = ecs.setPackingMass(packingMass);

            
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

            X = nan(4, 1);
            X(1 : 2) = ecs.getX(ne);
            X(3 : 4) = ecs.getX(pe);
            
        end

        function vals = getPhysicalValues(ecs, X)

            model = ecs.model;
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            
            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                
                vol  = sum(model.(elde).(am).G.cells.volumes);
                cmax = model.(elde).(am).(itf).cmax;

                vals.(elde).theta100       = X(2*ielde - 1);
                vals.(elde).alpha          = X(2*ielde)*vol*cmax;
                vals.(elde).volumeFraction = X(2*ielde);
                
            end
            
        end
        
        function X = getX(ecs, elde)

        % We define some shorthand names for simplicity.

            model = ecs.model;
            
            am  = 'ActiveMaterial';
            itf = 'Interface';

            X = nan(2, 1);
            
            X(1) = model.(elde).(am).(itf).theta100;
            vf   = unique(model.(elde).(am).volumeFraction)*model.(elde).(am).activeMaterialFraction;
            X(2) = vf;

        end
        
        
        function [fexp, fcomp] = setupfunction(ecs)

            fcomp = @(t, X) ecs.computeF(t, X);
            fexp  = @(t) ecs.experimentalF(t);
            
        end


        function c = conc(ecs, t, elde, theta, alpha)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';

            I = ecs.I;
            F = ecs.F;

            switch elde
              case ne
                sgn = -1;
              case pe
                sgn = 1;
            end
            
            c = theta + t * ((sgn*I)./(F*alpha));
            
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
            
            I     = ecs.I;
            model = ecs.model;
            T     = ecs.T;

            vals = ecs.getPhysicalValues(X);
            
            theta = vals.(pe).theta100;
            alpha = vals.(pe).alpha;

            c = ecs.conc(t, pe, theta, alpha);
            f = model.(pe).(am).(itf).computeOCPFunc(c, T, 1);
            
            theta = vals.(ne).theta100;
            alpha = vals.(ne).alpha;
            
            c = ecs.conc(t, ne, theta, alpha);
            f = f - model.(ne).(am).(itf).computeOCPFunc(c, T, 1);

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
            T         = ecs.T;
            totalTime = ecs.totalTime;
            
            props = ecs.getPhysicalValues(X);
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};
                
                cmax     = model.(elde).(am).(itf).cmax;
                theta100 = props.(elde).theta100;
                alpha    = props.(elde).alpha;
                
                theta0 = ecs.conc(totalTime, elde, theta100, alpha);
                
                switch elde
                  case ne
                    cM = theta100*cmax;
                    cm = theta0*cmax;
                  case pe
                    cM = theta0*cmax;
                    cm = theta100*cmax;
                end

                cap = (cM - cm)/cmax*alpha*F;
                
                props.(elde).theta0 = theta0;
                props.(elde).cmax   = cmax;
                props.(elde).cap    = cap;
                
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

            props.specEnergy = props.energy/ecs.mass;

            % props.NPratio = props.(ne).cap/props.(pe).cap;

            props.U = props.(pe).U - props.(ne).U;

        end

        function printParameters(ecs, X)

            vals = ecs.getPhysicalValues(X);
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';

            eldes = {ne, pe};

            fprintf('%20s%20s%20s\n', '', 'theta100', 'volume fraction');
            fprintf('%20s%20.5f%20.5f \n', ne, vals.(ne).theta100, vals.(ne).volumeFraction);
            fprintf('%20s%20.5f%20.5f \n', pe, vals.(pe).theta100, vals.(pe).volumeFraction);
            
        end

        function printProperties(ecs, X)

            props = ecs.computeProperties(X);

            packingMass = ecs.packingMass;
            mass        = ecs.mass;

            fprintf('%20s: %g [V]\n','initial voltage', props.U);
            fprintf('%20s: %g [Ah]\n', 'Capacity', props.cap/(1*hour));
            fprintf('%20s: %g [Wh/kg]\n', 'Specific Energy', props.specEnergy/(1*hour));
            % fprintf('%20s: %g [-]\n', 'N/P ratio', props.NPratio);
        end
        
    end


end
