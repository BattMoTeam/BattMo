classdef ParameterFitting
%% The goal of this class is to solve the following calibration problem:
% Given target values (for example NP ratio, total energy, see the complete list below),
% find the model parameters :
% - initial concentrations in electrodes (will be converted in theta100 values)
% - final concentrations in electrodes (will be converted in theta100 values)
% - active material volume fractions in electrodes
% - volume fractions in electrodes
    
    properties

        model
        con
        
        objective
        gradient
        constraints
        jacobian
        jacobianstructure

        lb % lower bounds on variables
        ub % upper bounds on variables
        cl % lower bounds on the constraints
        cu % upper bounds on the constraints

        % helper structures
        A

        cmaxs    % $c_{i,\max}$, maximum concentration
        Vs       % $V_i$, total volume of electrode (including pore space)

        packingMass
        mass
        masses
        
        targets % structure with field
                % - resLi       : residual factor for total amount of Lithium (see definition in method evaluateTargetValues)
                % - U           : Initial voltage
                % - energy      : battery energy
                % - specEnergy  : battery specific energy
                % - cap         : battery capacity
                % - specCap_neg : negative electrode specific capacity
                % - specCap_peg : positive electrode specific capacity
                % - vf_neg      : negative electrode volume fraction
                % - vf_pos      : positive electrode volume fraction
                % - NPratio     : Balancing coefficient
        
        weights % one value for each of the targets above

        used_targets % string list of target names 

        nparams % number of parameters
        
        ipopt % options that can be passed to ipopt
    end
    
    methods
        
        function pf = ParameterFitting(model, varargin)

            opt = struct('packingMass', 0);
            opt = merge_options(opt, varargin{:});
            
            pf.con   = PhysicalConstants();
            pf.model = model;
            
            pf = pf.setPackingMass(opt.packingMass);
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            Vs.n = sum(model.(ne).(co).G.getVolumes());
            Vs.p = sum(model.(pe).(co).G.getVolumes());

            cmaxs.n = model.(ne).(co).(am).(itf).saturationConcentration;
            cmaxs.p = model.(pe).(co).(am).(itf).saturationConcentration;

            % Targets with some example values
            targets.resLi      = 1;
            targets.U          = 4.2;
            targets.specEnergy = 280*watt*hour/kilogram;
            targets.cap        = 26.136*watt*hour;
            targets.NPratio    = 1.1;
            targets.avf_neg    = 0.9;
            targets.avf_pos    = 0.9;
            targets.vf_neg     = 0.75;
            targets.vf_pos     = 0.665;

            % Weight used for each of the target above
            weights.resLi      = 1;
            weights.U          = 1;
            weights.specEnergy = 1;
            weights.cap        = 1;
            weights.NPratio    = 1;
            weights.avf_neg    = 1;
            weights.avf_pos    = 1;
            weights.vf_neg     = 1;
            weights.vf_pos     = 1;
            
            weights = pf.normalizeWeights(weights);
            
            pf.Vs          = Vs;
            pf.cmaxs       = cmaxs;
            pf.targets     = targets;
            pf.weights     = weights;

            pf.used_targets = {'resLi', 'U', 'specEnergy', 'cap', 'NPratio', 'avf_neg', 'avf_pos', 'vf_neg', 'vf_pos'};
            pf.nparams = 8; % Number of parameters (see assignToX for overview)
            
            pf = pf.setupOptim();
            
        end

        function pf = setPackingMass(pf, packingMass)

            model = pf.model;
            
            [mass, masses] = computeCellMass(model, 'packingMass', packingMass);

            pf.packingMass = packingMass;
            pf.mass        = mass;
            pf.masses      = masses;
            
        end
        
        function weights = normalizeWeights(pf, weights)

            fdnames = pf.used_targets;
            
            tot = 0;
            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                tot = tot + weights.(fdname);
            end
            
            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                weights.(fdname) = 1/tot*weights.(fdname);
            end
            
        end
        
        function x = assignToX(pf, cs, avfs, vfs)
            
        % Variables
        % x(1) : $c_n^0$, initial concentration negative electrode
        % x(2) : $c_p^0$, initial concentration positive electrode
        % x(3) : $c_n^T$, concentration at final (discharge) time $T$ negative electrode
        % x(4) : $c_p^T$, concentration at final (discharge) time $T$ positive electrode
        % x(5) : $\gamma_i$, active material volume fraction negative electrode
        % x(6) : $\gamma_i$, active material volume fraction positiv electrode
        % x(7) : volume fraction negative electrode
        % x(8) : volume fraction positivt electrode
            
            x = nan(pf.nparams, 1);

            cmaxs = pf.cmaxs;
            
            x(1) = cs.n0/cmaxs.n;
            x(2) = cs.p0/cmaxs.p;
            x(3) = cs.nT/cmaxs.n;
            x(4) = cs.pT/cmaxs.p;
            x(5) = avfs.n;
            x(6) = avfs.p;
            x(7) = vfs.n;
            x(8) = vfs.p;
            
        end

        function [cs, avfs, vfs] = assignFromX(pf, x)

            if (nargin < 2) || isempty(x)
                x = pf.getModelValues();
            end

            cmaxs = pf.cmaxs;
            
            cs.n0  = x(1)*cmaxs.n;
            cs.p0  = x(2)*cmaxs.p;
            cs.nT  = x(3)*cmaxs.n;
            cs.pT  = x(4)*cmaxs.p;
            avfs.n = x(5);
            avfs.p = x(6);
            vfs.n  = x(7);
            vfs.p  = x(8);
            
        end

        function x = getModelValues(pf)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            model = pf.model;
            cmaxs = pf.cmaxs;
            
            avfs.n = model.(ne).(co).volumeFractions(1);
            vfs.n  = model.(ne).(co).volumeFraction;
            cs.n0  = model.(ne).(co).(am).(itf).guestStoichiometry100*cmaxs.n;
            cs.nT  = model.(ne).(co).(am).(itf).guestStoichiometry0*cmaxs.n;

            avfs.p = model.(pe).(co).volumeFractions(1);
            vfs.p  = model.(pe).(co).volumeFraction;
            cs.p0  = model.(pe).(co).(am).(itf).guestStoichiometry100*cmaxs.p;
            cs.pT  = model.(pe).(co).(am).(itf).guestStoichiometry0*cmaxs.p;

            x = pf.assignToX(cs, avfs, vfs);
            
        end
        
        function thetas = computeThetas(pf, x)
        % thetas.n0: $\theta_i^0$, initial lithiation
        % thetas.p0: $\theta_i^0$, initial lithiation
        % thetas.nT : $\theta_n^T$, lithiation at final (discharge) time $T$
        % thetas.pT : $\theta_p^T$, lithiation at final (discharge) time $T$

            if (nargin < 2) || isempty(x)
                x = pf.getModelValues();
            end

            cs    = pf.assignFromX(x);
            cmaxs = pf.cmaxs;

            thetas.n0 = cs.n0/cmaxs.n;
            thetas.p0 = cs.p0/cmaxs.p;
            thetas.nT = cs.nT/cmaxs.n;
            thetas.pT = cs.pT/cmaxs.p;            
            
        end

        function printParameters(pf, x)
        % Print the parameters, if x is not given or empty, use the parameters given by the model

            if (nargin < 2) || isempty(x)
                x = pf.getModelValues();
            end
            
            [cs, avfs, vfs] = pf.assignFromX(x);
            thetas          = pf.computeThetas(x);
            
            fprintf('%-60s: %g\n', 'theta start negative electrode', thetas.n0);
            fprintf('%-60s: %g\n', 'theta end negative electrode', thetas.nT);
            fprintf('%-60s: %g\n', 'theta start positive electrode', thetas.p0);
            fprintf('%-60s: %g\n', 'theta end positive electrode', thetas.pT);
            fprintf('%-60s: %g\n', 'active material volume fraction negative electrode', avfs.n);
            fprintf('%-60s: %g\n', 'active material volume fraction positive electrode', avfs.p);
            fprintf('%-60s: %g\n', 'volume fraction negative electrode', vfs.n);
            fprintf('%-60s: %g\n', 'volume fraction positive electrode', vfs.p);
            
        end

        function output = evaluateTargetValues(pf, x)
        % See output structure defined at the end of the function

            model    = pf.model;
            Vs       = pf.Vs;
            cmaxs    = pf.cmaxs;
            F        = pf.con.F;

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';

            if (nargin < 2) || isempty(x)
                x = pf.getModelValues();
            end

            [cs, avfs, vfs] = pf.assignFromX(x);

            eldes = {ne, pe};
            
            % setup helper vectors
            cM        = {cs.n0, cs.pT};
            cm        = {cs.nT, cs.p0};
            iVs       = {Vs.n, Vs.p};
            iavfs     = {avfs.n, avfs.p};
            ivolfracs = {vfs.n, vfs.p};
            
            for ielde = 1 : numel(eldes)
                
                n = 1; % charge number
                vol = iavfs{ielde}*ivolfracs{ielde}*iVs{ielde};
                caps{ielde} = (cM{ielde} - cm{ielde})*vol*n*F;
                
            end

            cap_neg = caps{1};
            cap_pos = caps{2};

            cap = cap_neg;
            ind = value(cap_neg) >= value(cap_pos);
            if any(ind)
                cap(ind) = cap_pos(ind);
            end

            c0s    = {cs.n0, cs.p0};
            cTs    = {cs.nT, cs.pT};
            icmaxs = {cmaxs.n, cmaxs.p};

            rs{1} = cap./cap_neg;
            rs{2} = cap./cap_pos;
            
            N = 1000;

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                
                smax = rs{ielde};
                s = smax.*linspace(0, 1, N + 1)';

                c = (1 - s).*c0s{ielde} + s.*cTs{ielde};
                f = model.(elde).(co).(am).(itf).computeOCP(c(1 : end - 1)/icmaxs{ielde});
                
                % function handler
                fs{ielde} = @(s) model.(elde).(co).(am).(itf).computeOCP(((1 - s).*c0s{ielde} + s.*cTs{ielde})/icmaxs{ielde});
                
                energies{ielde} = caps{ielde}*smax/N*sum(f);
                
            end
            
            energy = (energies{2} - energies{1});

            dischargeFunc = @(s) (fs{2}(s) - fs{1}(s));
            
            specEnergy = energy/pf.mass;

            specCap_neg = cap_neg/pf.masses.(ne).(co).val;
            specCap_pos = cap_pos/pf.masses.(pe).(co).val;

            nLi = Vs.n*avfs.n*vfs.n*cs.n0 + Vs.p*vfs.p*avfs.p*cs.p0;
            resLi = nLi./(Vs.p*vfs.p*avfs.p*cmaxs.p);

            NPratio = cap_neg./cap_pos;

            U = model.(pe).(co).(am).(itf).computeOCP(cs.p0, 298, cmaxs.p) - model.(ne).(co).(am).(itf).computeOCP(cs.n0/cmaxs.n);
            
            output = struct('nLi'          , nLi          , ...
                            'resLi'        , resLi        , ...
                            'U'            , U            , ...
                            'energy'       , energy       , ...
                            'dischargeFunc', dischargeFunc, ...
                            'specEnergy'   , specEnergy   , ...
                            'cap'          , cap          , ...
                            'cap_neg'      , cap_neg      , ...
                            'cap_pos'      , cap_pos      , ...
                            'specCap_neg'  , specCap_neg  , ...
                            'specCap_pos'  , specCap_pos  , ...
                            'avf_neg'      , avfs.n       , ...
                            'avf_pos'      , avfs.p       , ...
                            'vf_neg'       , vfs.n        , ...
                            'vf_pos'       , vfs.p        , ...
                            'NPratio'      , NPratio);
            
        end

        function printSpecifications(pf, x)
        % Print the targets, if x is not given or empty, use the parameters given by the model
            
            packingMass = pf.packingMass;
            mass        = pf.mass;

            if (nargin < 2) || isempty(x)
                x = pf.getModelValues();
            end
            
            output = pf.evaluateTargetValues(x);

            fprintf('%20s: %g [-]\n'    , 'N/P ratio'                 , output.NPratio);
            fprintf('%20s: %g [Ah]\n'   , 'Capacity'                  , output.cap/(1*hour));
            fprintf('%20s: %g [Wh/kg]\n', 'Specific Energy'           , output.specEnergy/(1*hour));
            fprintf('%20s: %g [mol]\n'  , 'total Lithium amount'      , output.nLi);
            fprintf('%20s: %g [V]\n'    , 'initial voltage'           , output.U);
            fprintf('%20s: %g [Ah]\n'   , 'Negative Capacity'         , output.cap_neg/(1*hour));
            fprintf('%20s: %g [Ah]\n'   , 'Positive Capacity'         , output.cap_pos/(1*hour));
            fprintf('%20s: %g [Ah/kg]\n', 'Negative Specific Capacity', output.specCap_neg/(1*hour));
            fprintf('%20s: %g [Ah/kg]\n', 'Positive Specific Capacity', output.specCap_pos/(1*hour));
            
        end


        function pf = setupOptim(pf)
            
        % Variables
        % x(1) : $c_n^0$, initial concentration negative electrode
        % x(2) : $c_p^0$, initial concentration positive electrode
        % x(3) : $c_n^T$, concentration at final (discharge) time $T$ negative electrode
        % x(4) : $c_p^T$, concentration at final (discharge) time $T$ positive electrode
        % x(5) : $\gamma_i$, active material volume fraction negative electrode
        % x(6) : $\gamma_i$, active material volume fraction positivt electrode
        % x(7) : volume fraction negative electrode
        % x(8) : volume fraction positivt electrode
            
        % Targets
        % - resLi       : residual factor for total amount of Lithium (see definition in method evaluateTargetValues)
        % - energy      : battery energy
        % - specEnergy  : battery specific energy
        % - U           : potential difference
        % - cap         : battery capacity
        % - specCap_neg : negative electrode specific capacity
        % - specCap_peg : positive electrode specific capacity
        % - vf_neg      : negative electrode specific capacity
        % - vf_peg      : positive electrode specific capacity
        % - NPratio     : Balancing coefficient

            cmaxs = pf.cmaxs;

            np = pf.nparams;
            
            lb = zeros(np, 1);
            ub = ones(np, 1);

            A = zeros(2, np);
            A(:, 1 : 4) = [[1 0 -1 0];
                           [0 -1 0 1]];
            A = sparse(A);
            cl = [0; 0];
            cu = [inf; inf];
            
            pf.lb = lb;
            pf.ub = ub;
            pf.cl = cl;
            pf.cu = cu;

            pf.A = A;

        end
        
        function y = objective_func(pf, x)
            
            o = pf.evaluateTargetValues(x);
            t = pf.targets;
            w = pf.weights;

            y = 0*x(1); % to get right AD initialization

            fdnames = pf.used_targets;
            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                y = y + w.(fdname)*(1/t.(fdname)*(o.(fdname) -t.(fdname)))^2;
            end

            
        end
        
        function y = gradient_func(pf, x)

            x = initVariablesADI(x);
            y = pf.objective_func(x);
            y = y.jac{1};
            y = y';
            
        end
        

        function y = constraints_func(pf, x)
            
            A = pf.A;
            y = A*x;

        end

        function y = jacobian_func(pf, x)

            x = initVariablesADI(x);
            y = pf.constraints_func(x);
            y = y.jac{1};
            
        end

        function y = jacobianstructure_func(pf)
        % we do not bother here and setup as full
            y = sparse(ones(2, pf.nparams));
        end


        function [x, info] = runOptim(pf)

            x0 = pf.getModelValues();

            funcs.objective         = @(x) pf.objective_func(x);
            funcs.gradient          = @(x) pf.gradient_func(x);
            funcs.constraints       = @(x) pf.constraints_func(x);            
            funcs.jacobian          = @(x) pf.jacobian_func(x);
            funcs.jacobianstructure = @() pf.jacobianstructure_func();
            
            options.lb = pf.lb;
            options.ub = pf.ub;
            options.cl = pf.cl;
            options.cu = pf.cu;
            
            options.ipopt = pf.ipopt;
            options.ipopt.hessian_approximation = 'limited-memory';

            % funcs.objective(x0)
            % funcs.gradient(x0)
            % funcs.constraints(x0)
            % funcs.jacobian(x0)
            % funcs.jacobianstructure()
            % return

            [x, info] = ipopt(x0, funcs, options);
            
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
