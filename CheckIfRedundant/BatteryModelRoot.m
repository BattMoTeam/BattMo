classdef BatteryModelRoot < RootModel

    properties

        con = physicalConstants();

        % Coupling structures
        couplingnames
        couplingTerms

        % fv2d structure (will disappear when we switch to own Newton solver)
        fv

        % Temperature and SOC
        % for the moment here, for convenience. Will be moved
        T
        SOC

        % Input current
        J
        % Voltage cut
        Ucut

    end

    methods

        function model = BatteryModelRoot(varargin)

            model = model@RootModel('battery');

            sepnx  = 30;
            nenx   = 30;
            penx   = 30;
            ccnenx = 20;
            ccpenx = 20;

            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            ny = 10;

            xlength = 1e-6*[10; 100; 50; 80; 10];
            ylength = 1e-2;

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            y = ylength/ny;
            y = rldecode(y, ny);
            y = [0; cumsum(y)];

            G = tensorGrid(x, y);
            G = computeGeometry(G);
            model.G = G;

            %% setup elyte
            nx = sum(nxs);

            submodels = {};

            istart = ccnenx + 1;
            ni = nenx + sepnx + penx;
            cells = pickTensorCells(istart, ni, nx, ny);
            submodels{end + 1} = orgLiPF6('elyte', G, cells);

            %% setup ne
            istart = ccnenx + 1;
            cells = pickTensorCells(istart, nenx, nx, ny);
            submodels{end + 1} = graphiteElectrode('ne', G, cells);

            %% setup pe
            istart = ccnenx + nenx + sepnx + 1;
            cells = pickTensorCells(istart, penx, nx, ny);
            submodels{end + 1} = nmc111Electrode('pe', G, cells);

            %% setup ccne
            istart = 1;
            cells = pickTensorCells(istart, ccnenx, nx, ny);
            submodels{end + 1} = currentCollector('ccne', G, cells);

            %% setup ccpe
            istart = ccnenx + nenx + sepnx + penx + 1;
            cells = pickTensorCells(istart, ccpenx, nx, ny);
            submodels{end + 1}  = currentCollector('ccpe', G, cells);

            %% setup sep
            istart = ccnenx + nenx + 1;
            cells = pickTensorCells(istart, sepnx, nx, ny);
            submodels{end + 1} = celgard2500('sep', G, cells);

            model.SubModels = submodels;
            model.hasparent = false;

            model = model.initiateCompositeModel();

            ccpe = model.getAssocModel('ccpe');
            ccpe.pnames = {ccpe.pnames{:}, 'E'};
            ccpe.names = {ccpe.names{:}, 'E'};
            model = model.setSubModel('ccpe', ccpe);

            model.names = {'T', 'SOC'};

            % setup ne

            ne = model.getAssocModel('ne');

            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            ne = ne.setPropFunction(PropFunction('T', fnupdate, fnmodel));
            ne = ne.setPropFunction(PropFunction('SOC', fnupdate, fnmodel));

            fnupdate = @(model, state) model.updatePhiElyte(state);
            fnmodel = {'..'};
            ne = ne.setPropFunction(PropFunction('phielyte', fnupdate, fnmodel));
            
            fnupdate = @(model, state) setupBCSources(model, state);
            fnmodel = {'..'};
            ne = ne.setPropFunction(PropFunction('jBcSource', fnupdate, fnmodel));

            model = model.setSubModel('ne', ne);

            % setup pe

            pe = model.getAssocModel('pe');

            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            pe = pe.setPropFunction(PropFunction('T', fnupdate, fnmodel));
            pe = pe.setPropFunction(PropFunction('SOC', fnupdate, fnmodel));

            fnupdate = @(model, state) model.updatePhiElyte(state);
            fnmodel = {'..'};
            pe = pe.setPropFunction(PropFunction('phielyte', fnupdate, fnmodel));

            fnupdate = @(model, state) setupBCSources(model, state);
            fnmodel = {'..'};
            pe = pe.setPropFunction(PropFunction('jBcSource', fnupdate, fnmodel));
            
            model = model.setSubModel('pe', pe);

            % setup ccne
            ccne = model.getAssocModel('ccne');
            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            ccne = ccne.setPropFunction(PropFunction('T', fnupdate, fnmodel));
            
            fnupdate = @(model, state) setupBCSources(model, state);
            fnmodel = {'..'};
            ccne = ccne.setPropFunction(PropFunction('jBcSource', fnupdate, fnmodel));
            
            model = model.setSubModel('ccne', ccne);

            % setup ccpe
            ccpe = model.getAssocModel('ccpe');
            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            ccpe = ccpe.setPropFunction(PropFunction('T', fnupdate, fnmodel));

            fnupdate = @(model, state) setupBCSources(model, state);
            fnmodel = {'..'};
            ccpe = ccpe.setPropFunction(PropFunction('jBcSource', fnupdate, fnmodel));
            
            model = model.setSubModel('ccpe', ccpe);

            % setup elyte
            elyte = model.getAssocModel('elyte');
            fnupdate = @(model, state) model.dispatchValues(state);
            fnmodel = {'..'};
            elyte = elyte.setPropFunction(PropFunction('T', fnupdate, fnmodel));

            model = model.setSubModel('elyte', elyte);

            model = model.initiateCompositeModel();

            model = model.setElytePorosity();

            %% setup couplings
            coupTerms = {};

            % coupling term 'ne-cc'
            coupTerms{end + 1} = setupNeElyteCoupTerm(model);
            coupTerms{end + 1} = setupPeElyteCoupTerm(model);
            coupTerms{end + 1} = setupCcneNeCoupTerm(model);
            coupTerms{end + 1} = setupCcpePeCoupTerm(model);
            coupTerms{end + 1} = setupCcneBcCoupTerm(model);
            coupTerms{end + 1} = setupCcpeBcCoupTerm(model);

            model.couplingTerms = coupTerms;
            model.couplingnames = cellfun(@(x) x.name, coupTerms, 'uniformoutput', false);

            model.SOC = 0.5;
            model.T   = 298.15;

            model.J = 0.1;
            model.Ucut = 2;

            model = model.initiateRootModel();
            
        end

        function state = initializeState(model, state)
        % Used only in debugging for the moment

            state = model.initiateState(state);

            % initialize each of the submodels
            for i = 1 : model.nSubModels
                state = model.SubModels{i}.initializeState(state);
            end

            ccne = model.getAssocModel('ccne');
            ne = model.getAssocModel('ne');

            [OCP_ne, state] = ne.getUpdatedProp(state, {'am', 'OCP'});
            nc = ccne.G.cells.num;
            OCP_ccne = OCP_ne(1)*ones(nc, 1);

            state = ccne.setProp(state, 'OCP', OCP_ccne);
            state = ccne.setProp(state, 'phi', OCP_ccne);

            ne = model.getAssocModel('ccpe');
            pe = model.getAssocModel('pe');
            ccpe = model.getAssocModel('ccpe');

            [OCP_pe, state] = pe.getUpdatedProp(state, {'am', 'OCP'});
            nc = ccpe.G.cells.num;
            OCP_ccpe = OCP_pe(1)*ones(nc, 1);

            state = ccpe.setProp(state, 'OCP', OCP_ccpe);
            state = ccpe.setProp(state, 'phi', OCP_ccpe);

            state = ccpe.setProp(state, 'E', OCP_ccpe(1));

        end

        function state = dispatchValues(model, state)

            [T, state] = model.getUpdatedProp(state, 'T');
            [SOC, state] = model.getUpdatedProp(state, 'SOC');

            elyte = model.getAssocModel('elyte');
            G = elyte.G;
            Telyte = T(G.mappings.cellmap);

            ne = model.getAssocModel('ne');
            G = ne.G;
            Tne = T(G.mappings.cellmap);
            SOCne = SOC(G.mappings.cellmap);

            pe = model.getAssocModel('pe');
            G = pe.G;
            Tpe = T(G.mappings.cellmap);
            SOCpe = SOC(G.mappings.cellmap);

            ccpe = model.getAssocModel('ccpe');
            G = ccpe.G;
            Tccpe = T(G.mappings.cellmap);

            ccne = model.getAssocModel('ccne');
            G = ccne.G;
            Tccne = T(G.mappings.cellmap);

            state = model.setProp(state, {'elyte', 'T'}, Telyte);
            state = model.setProp(state, {'ne', 'T'}, Tne);
            state = model.setProp(state, {'pe', 'T'}, Tpe);
            state = model.setProp(state, {'ccpe', 'T'}, Tccpe);
            state = model.setProp(state, {'ccne', 'T'}, Tccne);
            state = model.setProp(state, {'ne', 'SOC'}, SOCne);
            state = model.setProp(state, {'pe', 'SOC'}, SOCpe);

        end

        function state = updatePhiElyte(model, state)

            [phielyte, state] = model.getUpdatedProp(state, {'elyte', 'phi'});
            elyte = model.getAssocModel('elyte');

            elytecells = zeros(model.G.cells.num, 1);
            elytecells(elyte.G.mappings.cellmap) = (1 : elyte.G.cells.num)';

            ne = model.getAssocModel('ne');
            phielyte_ne = phielyte(elytecells(ne.G.mappings.cellmap));

            pe = model.getAssocModel('pe');
            phielyte_pe = phielyte(elytecells(pe.G.mappings.cellmap));

            state = model.setProp(state, {'ne', 'phielyte'}, phielyte_ne);
            state = model.setProp(state, {'pe', 'phielyte'}, phielyte_pe);

        end

        
        function model = setElytePorosity(model)

            elyte = model.getAssocModel('elyte');
            ne = model.getAssocModel('ne');
            pe = model.getAssocModel('pe');
            sep = model.getAssocModel('sep');

            elytecells = zeros(model.G.cells.num, 1);
            elytecells(elyte.G.mappings.cellmap) = (1 : elyte.G.cells.num)';

            elyte.eps = NaN(elyte.G.cells.num, 1);
            elyte.eps(elytecells(ne.G.mappings.cellmap)) = ne.void;
            elyte.eps(elytecells(pe.G.mappings.cellmap)) = pe.void;
            elyte.eps(elytecells(sep.G.mappings.cellmap)) = sep.void;

            model = model.setSubModel('elyte', elyte);

        end

        function model = setupFV(model, state)
            model.fv = fv2d(model, state);
        end

        function [y0, yp0] = dynamicPreprocess(model, state)
            %% Initialize the state vector

            y0 = [];
            pvarnames = model.getModelPrimaryVarNames();
            for ind = 1 : numel(pvarnames);
                y0 = [y0; model.getProp(state, pvarnames{ind})];
            end

            yp0  = zeros(length(y0), 1);

        end

        function [t, y] = p2d(model)

            % Set initial conditions
            initstate = model.icp2d();

            % Generate the FV structure
            model.fv = fv2d(model, initstate);

            %% Solution space (time or current) discretization

            % Time discretization
            model.fv.ti = 0;
            model.fv.tf = 3600*24;
            model.fv.dt = 10;
            model.fv.tUp = 0.1;
            model.fv.tSpan = (model.fv.ti:model.fv.dt:model.fv.tf);

            % Pre-process
            [y0, yp0] = model.dynamicPreprocess(initstate);

            % Set up and solve the system of equations
            endFun = @(t,y,yp) model.cutOff(t, y, yp);
            fun    = @(t,y,yp) model.odefun(t, y, yp);
            derfun = @(t,y,yp) model.odederfun(t, y, yp);

            options = odeset('RelTol'  , 1e-4  , ...
                             'AbsTol'  , 1e-6  , ...
                             'Stats'   , 'on'  , ...
                             'Events'  , endFun, ...
                             'Jacobian', derfun);

            [t, y] = ode15i(fun, model.fv.tSpan', y0, yp0, options);

        end

        function res = odefun(model, t, y, yp, varargin)
        %ODEFUN Compiles the system of differential equations

            opt = struct('useAD', false);
            opt = merge_options(opt, varargin{:});
            useAD = opt.useAD;

            % Build SOE
            res = model.dynamicBuildSOE(t, y, yp, 'useAD', useAD);
            
        end

        function [dfdy, dfdyp] = odederfun(model, t, y, yp)
            
            res = model.odefun(t, y, yp, 'useAD', true);
            if(numel(res.jac) == 2)
                dfdy  = res.jac{1};
                dfdyp = res.jac{2};
            else
                dfdy  = res.jac{1}(:,1:res.numVars(1));
                dfdyp = res.jac{1}(:,res.numVars(1)+1:end);
            end
            
        end

        function [value, isterminal, direction] = cutOff(model, t, y, yp)
        % This will be reimplemented when we move away from ode15i

            varname = VarName({'ccpe'}, 'E');

            % here we assume E at ccne is equal to zero
            U = y(model.fv.getSlot(varname));

            value = U - model.Ucut;
            isterminal = 1;
            direction = 0;
        end

        function initstate = icp2d(model)
        % Setup initial state

            nc = model.G.cells.num;

            initstate = [];
            initstate = model.initiateState(initstate);

            SOC = model.SOC;
            T   = model.T;
            initstate = model.setProp(initstate, 'T', T*ones(nc, 1));
            initstate = model.setProp(initstate, 'SOC', SOC*ones(nc, 1));

            elyte = model.getAssocModel('elyte');
            ne    = model.getAssocModel('ne');
            pe    = model.getAssocModel('pe');
            ccne  = model.getAssocModel('ccne');
            ccpe  = model.getAssocModel('ccpe');

            ne_am = ne.getAssocModel('am');
            pe_am = pe.getAssocModel('am');

            %% setup initial ne state

            m = (1 ./ (ne_am.theta100 - ne_am.theta0));
            b = -m .* ne_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* ne_am.Li.cmax;
            c = c*ones(ne.G.cells.num, 1);

            initstate = ne_am.setProp(initstate, 'Li', c);
            [OCP, initstate] = ne_am.getUpdatedProp(initstate, 'OCP');
            initstate = ne_am.setProp(initstate, 'phi', OCP);

            %% setup initial pe state

            m = (1 ./ (pe_am.theta100 - pe_am.theta0));
            b = -m .* pe_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* pe_am.Li.cmax;
            c = c*ones(ne.G.cells.num, 1);

            initstate = pe_am.setProp(initstate, 'Li', c);
            [OCP, initstate] = pe_am.getUpdatedProp(initstate, 'OCP');
            initstate = pe_am.setProp(initstate, 'phi', OCP);

            %% setup initial elyte state

            initstate = elyte.setProp(initstate, 'phi', zeros(elyte.G.cells.num, 1));
            initstate = elyte.setProp(initstate, 'cLi', 1000*ones(elyte.G.cells.num, 1));

            %% setup initial Current collectors state
            [OCP, initstate] = ne_am.getUpdatedProp(initstate, 'OCP');
            OCP = OCP(1) .* ones(ccne.G.cells.num, 1);
            initstate = ccne.setProp(initstate, 'phi', OCP);

            [OCP, initstate] = pe_am.getUpdatedProp(initstate, 'OCP');
            OCP = OCP(1) .* ones(ccpe.G.cells.num, 1);
            initstate = ccpe.setProp(initstate, 'phi', OCP);
            initstate = ccpe.setProp(initstate, 'E', OCP(1));
            
        end

        function soe = dynamicBuildSOE(model, t, y, yp, varargin)

            opt = struct('useAD', false);
            opt = merge_options(opt, varargin{:});
            useAD = opt.useAD;

            fv = model.fv;

            if useAD
                adbackend = model.AutoDiffBackend();
                [y, yp] = adbackend.initVariablesAD(y, yp);
            end

            % Mapping of variables
            state = [];
            state = model.initiateState(state);
            nc = model.G.cells.num;

            % setup temperature and SOC here
            SOC = model.SOC;
            T   = model.T;
            state = model.setProp(state, 'T', T*ones(nc, 1));
            state = model.setProp(state, 'SOC', SOC*ones(nc, 1));

            pvarnames = model.getModelPrimaryVarNames();
            for ind = 1 : numel(pvarnames)
                varname = pvarnames{ind};
                state = model.setProp(state, varname, y(fv.getSlot(varname)));
            end

            % variables for time derivatives
            elyte_Li_cdot = yp(fv.getSlot(VarName({'elyte'}, 'cLi')));
            ne_Li_csdot   = yp(fv.getSlot(VarName({'ne', 'am'}, 'Li')));
            pe_Li_csdot   = yp(fv.getSlot(VarName({'pe', 'am'}, 'Li')));

            elyte = model.getAssocModel('elyte');
            ne    = model.getAssocModel('ne');
            pe    = model.getAssocModel('pe');
            ccne  = model.getAssocModel('ccne');
            ccpe  = model.getAssocModel('ccpe');

            ne_am = ne.getAssocModel('am');
            pe_am = pe.getAssocModel('am');

            elyte_cLi = elyte.getUpdatedProp(state, 'cLi');
            elyte_phi  = elyte.getUpdatedProp(state, 'phi');
            ne_Li      = ne.getUpdatedProp(state, {'am', 'Li'});
            ne_phi     = ne.getUpdatedProp(state, {'am', 'phi'});
            pe_Li      = pe.getUpdatedProp(state, {'am', 'Li'});
            pe_phi     = pe.getUpdatedProp(state, {'am', 'phi'});
            ccne_phi   = ccne.getUpdatedProp(state, 'phi');
            ccpe_phi   = ccpe.getUpdatedProp(state, 'phi');
            ccpe_E     = ccpe.getUpdatedProp(state, 'E');

            elyte_Li_cepsdot = elyte.eps.*elyte_Li_cdot;
            ne_Li_csepsdot   = ne_am.eps.*ne_Li_csdot;
            pe_Li_csepsdot   = pe_am.eps.*pe_Li_csdot;



            %% Cell voltage
            ccne_E = 0;
            ccpe_E = ccpe.getProp(state, 'E');

            U = ccpe_E - ccne_E;

            if useAD
                adsample = getSampleAD(y, yp);
                adbackend = model.AutoDiffBackend;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Source terms for continuity equations                    %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % We setup the exchange terms between the electrolyte and the electrodes for Li due to chemical reacions:
            % elyte_Li_source, ne_Li_source, pe_Li_source.
            %
            % We setup also the electron production in the electrod: ne_e_source, pe_e_source.
            %

            % Inititalize vectors
            elyte_Li_source = zeros(elyte.G.cells.num, 1);
            ne_Li_source    = zeros(ne.G.cells.num, 1);
            pe_Li_source    = zeros(pe.G.cells.num, 1);
            ne_e_source     = zeros(ne.G.cells.num, 1);
            pe_e_source     = zeros(pe.G.cells.num, 1);

            if useAD
                elyte_Li_source = adbackend.convertToAD(elyte_Li_source, adsample);
                ne_Li_source    = adbackend.convertToAD(ne_Li_source, adsample);
                pe_Li_source    = adbackend.convertToAD(pe_Li_source, adsample);
                ne_e_source     = adbackend.convertToAD(ne_e_source, adsample);
                pe_e_source     = adbackend.convertToAD(pe_e_source, adsample);
            end

            %%%%% Set up chemical source terms %%%%%%%%%%%%%%%%%k

            %%%%% NE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            coupterm = model.getCoupTerm('ne-elyte');
            necells = coupterm.couplingcells(:, 1);
            elytecells = coupterm.couplingcells(:, 2);

            % calculate rection rate
            [ne_R, state] = ne.getUpdatedProp(state, 'R');

            % Electrolyte NE Li+ source
            elyte_Li_source(elytecells) = ne_R;

            % Active Material NE Li0 source
            ne_Li_source(necells) = - ne_R;

            % Active Material NE current source
            ne_e_source(necells) = + ne_R;

            %%%%% PE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Electrolyte PE Li+ source
            coupterm = model.getCoupTerm('pe-elyte');
            pecells = coupterm.couplingcells(:, 1);
            elytecells = coupterm.couplingcells(:, 2);

            % calculate rection rate
            [pe_R, state] = pe.getUpdatedProp(state, 'R');

            % Electrolyte PE Li+ source
            elyte_Li_source(elytecells) = - pe_R;

            % Active Material PE Li0 source
            pe_Li_source(pecells) = + pe_R;

            % Active Material PE current source
            pe_e_source(pecells) = - pe_R;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Diffusion Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of diffusion mass flux

            [D, state] = elyte.getUpdatedProp(state, 'D');
            elyte_Deff = D .* elyte.eps .^1.5;


            x = elyte.getUpdatedProp(state, 'cLi');
            trans = elyte.operators.harmFace(elyte_Deff);
            flux = - trans.*elyte.operators.Grad(x);
            elyte_Li_divDiff = elyte.operators.Div(flux)./elyte.G.cells.volumes;

            x = ne.getUpdatedProp(state, {'am', 'Li'});
            trans = ne.operators.harmFace(ne_Deff);
            flux = - trans.*ne.operators.Grad(x);
            ne_Li_divDiff = ne.operators.Div(flux)./ne.G.cells.volumes;

            x = pe.getUpdatedProp(state, {'am', 'Li'});
            trans = pe.operators.harmFace(pe_Deff);
            flux = - trans.*pe.operators.Grad(x);
            pe_Li_divDiff = pe.operators.Div(flux)./pe.G.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Migration Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of the migration mass flux for the electrolyte Li+ Migration

            [elyte_j, state] = elyte.getUpdatedProp(state, 'j');
            ind_Li = 1;
            flux = elyte.sp.t{ind_Li} ./ (elyte.sp.z{ind_Li} .* model.con.F) .* elyte_j;
            elyte_Li_divMig = elyte.operators.Div(flux)./elyte.G.cells.volumes;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% System of Equations                                      %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Liquid electrolyte dissolved ionic species mass continuity %
            %   Electrolyte Li+ Mass Continuity
            elyte_Li_massCont = (-elyte_Li_divDiff - elyte_Li_divMig + elyte_Li_source - elyte_Li_cepsdot);

            %% Liquid electrolyte charge continuity %%%%%%%%%%%%%%%%%%%%%%%

            elyte_chargeCons = -(elyte.operators.Div(elyte_j)./elyte.G.cells.volumes)./model.con.F + ...
                elyte_Li_source.*elyte.sp.z{ind_Li};

            %% Electrode Active material mass continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%

            ne_Li_massCont = (-ne_Li_divDiff + ne_Li_source - ne_Li_csepsdot);
            pe_Li_massCont = (-pe_Li_divDiff + pe_Li_source - pe_Li_csepsdot);

            %% Electode Active material charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%

            [ne_j, state] = ne.getUpdatedProp(state, 'j');
            [pe_j, state] = pe.getUpdatedProp(state, 'j');
            [ne_j_bcsource, state] = ne.getUpdatedProp(state, 'jBcSource');
            [pe_j_bcsource, state] = pe.getUpdatedProp(state, 'jBcSource');
            
            ne_e_chargeCons = (ne.operators.Div(ne_j) - ne_j_bcsource)./ ne.G.cells.volumes./model.con.F - ...
                ne_e_source;
            pe_e_chargeCons = (pe.operators.Div(pe_j) - pe_j_bcsource)./ pe.G.cells.volumes./model.con.F - ...
                pe_e_source;

            %% Collector charge continuity

            [ccne_j, state] = ccne.getUpdatedProp(state, 'j');
            [ccpe_j, state] = ccpe.getUpdatedProp(state, 'j');
            [ccne_j_bcsource, state] = ccne.getUpdatedProp(state, 'jBcSource');
            [ccpe_j_bcsource, state] = ccpe.getUpdatedProp(state, 'jBcSource');
            
            ccne_e_chargeCons = (ccne.operators.Div(ccne_j) - ccne_j_bcsource)./ ccne.G.cells.volumes./model.con.F;
            ccpe_e_chargeCons = (ccpe.operators.Div(ccpe_j) - ccpe_j_bcsource)./ ccpe.G.cells.volumes./model.con.F;

            %% Control equation

            src = currentSource(t, fv.tUp, fv.tf, model.J);
            coupterm = model.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = ccpe_E;
            [tccpe, cells] = ccpe.operators.harmFaceBC(ccpe_sigmaeff, faces);
            control = src - sum(tccpe.*(bcval - ccpe_phi(cells)));

            %% Governing equations

            soe = vertcat(elyte_Li_massCont, ...
                          elyte_chargeCons , ...
                          ne_Li_massCont   , ...
                          ne_e_chargeCons  , ...
                          pe_Li_massCont   , ...
                          pe_e_chargeCons  , ...
                          ccne_e_chargeCons, ...
                          ccpe_e_chargeCons, ...
                          control);

        end

        function coupterm = getCoupTerm(model, coupname)
            coupnames = model.couplingnames;

            [isok, ind] = ismember(coupname, coupnames);
            assert(isok, 'name of coupling term is not recognized.');

            coupterm = model.couplingTerms{ind};

        end

        function coupTerm = setupNeElyteCoupTerm(model)

            ne = model.getAssocModel('ne');
            elyte = model.getAssocModel('elyte');

            Gne = ne.G;
            Gelyte = elyte.G;

            % parent Grid
            G = Gne.mappings.parentGrid;

            % All the cells from ne are coupled with elyte
            cells1 = (1 : Gne.cells.num)';
            pcells = Gne.mappings.cellmap(cells1);

            mapping = zeros(G.cells.num, 1);
            mapping(Gelyte.mappings.cellmap) = (1 : Gelyte.cells.num)';
            cells2 = mapping(pcells);

            compnames = {'ne', 'elyte'};
            coupTerm = couplingTerm('ne-elyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty

        end

        function coupTerm = setupPeElyteCoupTerm(model)

            pe = model.getAssocModel('pe');
            elyte = model.getAssocModel('elyte');

            Gpe = pe.G;
            Gelyte = elyte.G;

            % parent Grid
            G = Gpe.mappings.parentGrid;

            % All the cells from pe are coupled with elyte
            cells1 = (1 : Gpe.cells.num)';
            pcells = Gpe.mappings.cellmap(cells1);

            mapping = zeros(G.cells.num, 1);
            mapping(Gelyte.mappings.cellmap) = (1 : Gelyte.cells.num)';
            cells2 = mapping(pcells);

            compnames = {'pe', 'elyte'};
            coupTerm = couplingTerm('pe-elyte', compnames);
            coupTerm.couplingcells = [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling between faces

        end

        function coupTerm = setupCcneNeCoupTerm(model)

            ne = model.getAssocModel('ne');
            ccne = model.getAssocModel('ccne');

            Gne = ne.G;
            Gccne = ccne.G;

            % parent Grid
            G = Gne.mappings.parentGrid;

            % We pick up the faces at the right of Cccne
            xf = Gccne.faces.centroids(:, 1);
            mxf = max(xf);
            faces1 = find(xf > (1 - eps)*mxf);

            pfaces = Gccne.mappings.facemap(faces1);
            mapping = zeros(G.faces.num, 1);
            mapping(Gne.mappings.facemap) = (1 : Gne.faces.num)';
            faces2 = mapping(pfaces);

            cells1 = sum(Gccne.faces.neighbors(faces1, :), 2);
            cells2 = sum(Gne.faces.neighbors(faces2, :), 2);

            compnames = {'ccne', 'ne'};
            coupTerm = couplingTerm('ccne-ne', compnames);
            coupTerm.couplingfaces =  [faces1, faces2];
            coupTerm.couplingcells = [cells1, cells2];

        end

        function coupTerm = setupCcpePeCoupTerm(model)

            pe = model.getAssocModel('pe');
            ccpe = model.getAssocModel('ccpe');

            Gpe = pe.G;
            Gccpe = ccpe.G;

            % parent Grid
            G = Gpe.mappings.parentGrid;

            % We pick up the faces at the left of Cccpe
            xf = Gccpe.faces.centroids(:, 1);
            mxf = min(xf);
            faces1 = find(xf < (1 + eps)*mxf);

            pfaces = Gccpe.mappings.facemap(faces1);
            mapping = zeros(G.faces.num, 1);
            mapping(Gpe.mappings.facemap) = (1 : Gpe.faces.num)';
            faces2 = mapping(pfaces);

            cells1 = sum(Gccpe.faces.neighbors(faces1, :), 2);
            cells2 = sum(Gpe.faces.neighbors(faces2, :), 2);

            compnames = {'ccpe', 'pe'};
            coupTerm = couplingTerm('ccpe-pe', compnames);
            coupTerm.couplingfaces =  [faces1, faces2];
            coupTerm.couplingcells = [cells1, cells2];

        end

        function coupTerm = setupCcneBcCoupTerm(model)

            ccne = model.getAssocModel('ccne');
            G = ccne.G;

            % We pick up the faces at the top of Cccne
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(yf > (1 - eps)*myf);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccne'};
            coupTerm = couplingTerm('bc-ccne', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

        function coupTerm = setupCcpeBcCoupTerm(model)

            ccpe = model.getAssocModel('ccpe');
            G = ccpe.G;

            % We pick up the faces at the top of Cccpe
            yf = G.faces.centroids(:, 2);
            myf = max(yf);
            faces = find(yf > (1 - eps)*myf);
            cells = sum(G.faces.neighbors(faces, :), 2);

            compnames = {'ccpe'};
            coupTerm = couplingTerm('bc-ccpe', compnames);
            coupTerm.couplingfaces = faces;
            coupTerm.couplingcells = cells;

        end

    end

end
