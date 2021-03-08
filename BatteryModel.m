classdef BatteryModel < CompositeModel

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

        % submodels
        elyte
        ne
        pe
        ccpe
        ccne
        
    end

    methods

        function model = BatteryModel(varargin)

            model = model@CompositeModel('battery');
            model.AutoDiffBackend =  SparseAutoDiffBackend('useBlocks',false);
            names = {'T', 'SOC'};
            model.names = names;
            model = model.setupVarDims();

            model = model.setupBatteryComponents();
            
            model.hasparent = false;

            model = model.initiateCompositeModel();

            ccpe = model.getAssocModel('ccpe');
            ccpe.pnames = {ccpe.pnames{:}, 'E'};
            ccpe.names = {ccpe.names{:}, 'E'};
            ccpe.vardims('E') = 1;
            model = model.setSubModel('ccpe', ccpe);

            model.names = {'T', 'SOC'};

            %% Setup property update functions 

            ne = model.getAssocModel('ne');
            pe = model.getAssocModel('pe');
            ccne = model.getAssocModel('ccne');
            ccpe = model.getAssocModel('ccpe');
            elyte = model.getAssocModel('elyte');
            
            % update function for temperature and soc
            fnupdate = @(model, state) model.dispatchValues(state);
            inputnames = {VarName({'..'}, 'T'), ...
                          VarName({'..'}, 'SOC')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('T'  , fnupdate, inputnames, fnmodel);
            ne = ne.addPropFunction('SOC', fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('T'  , fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('SOC', fnupdate, inputnames, fnmodel);
            ccne = ccne.addPropFunction('T', fnupdate, inputnames, fnmodel);
            ccpe = ccpe.addPropFunction('T', fnupdate, inputnames, fnmodel);
            elyte = elyte.addPropFunction('T', fnupdate, inputnames, fnmodel);
            
            % update function for exchange term (ne-elyte)
            fnupdate = @(model, state) setupExchanges(model, state);
            inputnames = {VarName({'..', 'ne'}, 'R'), ...
                          VarName({'..', 'pe'}, 'R')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('LiSource', fnupdate, inputnames, fnmodel);
            ne = ne.addPropFunction('eSource' , fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('LiSource', fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('eSource' , fnupdate, inputnames, fnmodel);
            elyte = elyte.addPropFunction('LiSource', fnupdate, inputnames, fnmodel);
            
            % update function for phielyte (electrolyte potential)
            fnupdate = @(model, state) model.updatePhiElyte(state);
            inputnames = {VarName({'..', 'elyte'}, 'phi')};
            fnmodel = {'..'};
            ne = ne.addPropFunction('phielyte', fnupdate, inputnames, fnmodel);
            pe = pe.addPropFunction('phielyte', fnupdate, inputnames, fnmodel);
            
            % update function for boundary terms (ne-ccne)
            fnupdate = @(model, state) setupBCSources(model, state);
            inputnames = {VarName({'..', 'ne', 'am'}, 'phi'), ...
                          VarName({'..', 'pe', 'am'}, 'phi'), ... 
                          VarName({'..', 'ccne'}, 'phi'), ...
                          VarName({'..', 'ccpe'}, 'phi'), ...
                          VarName({'..', 'ccpe'}, 'E'), ...
                         };
            fnmodel = {'..'};
            ne   = ne.addPropFunction('jBcSource', fnupdate, inputnames, fnmodel);
            pe   = pe.addPropFunction('jBcSource', fnupdate, inputnames, fnmodel);
            ccne = ccne.addPropFunction('jBcSource', fnupdate, inputnames, fnmodel);
            ccpe = ccpe.addPropFunction('jBcSource', fnupdate, inputnames, fnmodel);
            
            
            fnupdate = @(model, state) model.dynamicBuildSOE(state); 
            % function above is not correct. This is just used now to inform the graph
            inputnames = {VarName({'..', 'ne'}   , 'LiSource'), ...
                          VarName({'..', 'pe'}   , 'LiSource'), ...
                          VarName({'..', 'elyte'}, 'LiSource'), ...
                          VarName({'..', 'ne'}   , 'LiFlux'), ...
                          VarName({'..', 'pe'}   , 'LiFlux'), ...
                          VarName({'..', 'elyte'}, 'LiFlux'), ...
                         };
            fnmodel = {'..'};
            ne    = ne.addPropFunction('massCont', fnupdate, inputnames, fnmodel);
            pe    = pe.addPropFunction('massCont', fnupdate, inputnames, fnmodel);
            elyte = elyte.addPropFunction('massCont', fnupdate, inputnames, fnmodel);
            
            model = model.setSubModel('ne', ne);
            model = model.setSubModel('pe', pe);
            model = model.setSubModel('ccne', ccne);
            model = model.setSubModel('ccne', ccne);
            model = model.setSubModel('elyte', elyte);

            model = model.initiateCompositeModel();

            %% setup porosity in all components
            model = model.setElytePorosity();
            
            
            model = model.initiateCompositeModel();
            
            model.elyte = model.getAssocModel('elyte');
            model.ne    = model.getAssocModel('ne');
            model.pe    = model.getAssocModel('pe');
            model.ccne  = model.getAssocModel('ccne');
            model.ccpe  = model.getAssocModel('ccpe');
            
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
            model.T = 298.15;

            model.J = 0.1;
            model.Ucut = 2;
            
        end
        
        
        function model = setupBatteryComponents(model)
            
            fac = 1;
            sepnx  = 30*fac;
            nenx   = 30*fac;
            penx   = 30*fac;
            ccnenx = 20*fac;
            ccpenx = 20*fac;

            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            ny = 10*fac;

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
            
        end
        
        function state = dispatchValues(model, state)

            T = state.T;
            SOC = state.SOC;

            elyte = model.elyte;
            G = elyte.G;
            Telyte = T(G.mappings.cellmap);

            ne = model.ne;
            G = ne.G;
            Tne = T(G.mappings.cellmap);
            SOCne = SOC(G.mappings.cellmap);
            
            pe = model.pe;
            G = pe.G;
            Tpe = T(G.mappings.cellmap);
            SOCpe = SOC(G.mappings.cellmap);
            
            ccpe = model.ccpe;
            G = ccpe.G;
            Tccpe = T(G.mappings.cellmap);

            ccne = model.ccne;
            G = ccne.G;
            Tccne = T(G.mappings.cellmap);

            state.elyte.T = Telyte;
            state.ne.T    = Tne;
            state.ne.am.T = Tne;
            state.pe.T    = Tpe;
            state.pe.am.T = Tpe;
            state.ccpe.T  = Tccpe;
            state.ccne.T  = Tccne;
            
            state.ne.SOC    = SOCne;
            state.ne.am.SOC = SOCne;
            state.pe.SOC    = SOCpe;
            state.pe.am.SOC = SOCpe;
            
        end

        function state = updatePhiElyte(model, state)

            phielyte = state.elyte.phi;
            elyte = model.elyte;

            elytecells = zeros(model.G.cells.num, 1);
            elytecells(elyte.G.mappings.cellmap) = (1 : elyte.G.cells.num)';

            ne = model.ne;
            phielyte_ne = phielyte(elytecells(ne.G.mappings.cellmap));

            pe = model.pe;
            phielyte_pe = phielyte(elytecells(pe.G.mappings.cellmap));

            state.ne.phielyte = phielyte_ne;
            state.pe.phielyte = phielyte_pe;

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
            
            y0   = [ state.elyte.cs{1};
                     state.elyte.phi;
                     state.ne.am.Li;
                     state.ne.am.phi;
                     state.pe.am.Li;
                     state.pe.am.phi;
                     state.ccne.phi;
                     state.ccpe.phi;
                     state.ccpe.E];
            
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
            %model.fv.tSpan = (model.fv.ti:model.fv.dt:model.fv.tf);
            model.fv.tSpan = [model.fv.ti,model.fv.tf];

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

            % here we assume E at ccne is equal to zero
            U = y(model.fv.slots{end});

            value = U - model.Ucut;
            isterminal = 1;
            direction = 0;
        end

        function initstate = icp2d(model)
        % Setup initial state

            nc = model.G.cells.num;

            SOC = model.SOC;
            T   = model.T;
            
            initstate.T   =  T*ones(nc, 1);
            initstate.SOC =  SOC*ones(nc, 1);

            initstate = model.dispatchValues(initstate);
            
            elyte = model.elyte;
            ne    = model.ne;
            pe    = model.pe;
            ccne  = model.ccne;
            ccpe  = model.ccpe;

            ne_am = ne.am;
            pe_am = pe.am;

            %% setup initial ne state

            m = (1 ./ (ne_am.theta100 - ne_am.theta0));
            b = -m .* ne_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* ne_am.Li.cmax;
            c = c*ones(ne.G.cells.num, 1);

            initstate.ne.am.Li = c;
            initstate.ne.am = ne_am.updateQuantities(initstate.ne.am);

            OCP = initstate.ne.am.OCP;
            initstate.ne.am.phi = OCP;

            %% setup initial pe state

            m = (1 ./ (pe_am.theta100 - pe_am.theta0));
            b = -m .* pe_am.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* pe_am.Li.cmax;
            c = c*ones(ne.G.cells.num, 1);

            initstate.pe.am.Li = c;
            initstate.pe.am = pe_am.updateQuantities(initstate.pe.am);

            OCP = initstate.pe.am.OCP;
            initstate.pe.am.phi = OCP;

            %% setup initial elyte state

            initstate.elyte.phi = zeros(elyte.G.cells.num, 1);
            cs=cell(2,1);
            initstate.elyte.cs=cs;
            initstate.elyte.cs{1} = 1000*ones(elyte.G.cells.num, 1);

            %% setup initial Current collectors state
            
            OCP = initstate.ne.am.OCP;
            OCP = OCP(1) .* ones(ccne.G.cells.num, 1);
            initstate.ccne.phi = OCP;

            OCP = initstate.pe.am.OCP;
            OCP = OCP(1) .* ones(ccpe.G.cells.num, 1);
            initstate.ccpe.phi = OCP;
            
            initstate.ccpe.E = OCP(1);
            
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
            nc = model.G.cells.num;

            % setup temperature and SOC here
            SOC = model.SOC;
            T   = model.T;
            state.T =  T*ones(nc, 1);
            state.SOC =  SOC*ones(nc, 1);

            sl = fv.slots;
            
            state.elyte.cs{1} = y(sl{1});
            state.elyte.phi   = y(sl{2});
            state.ne.am.Li    = y(sl{3});
            state.ne.am.phi   = y(sl{4});
            state.pe.am.Li    = y(sl{5});
            state.pe.am.phi   = y(sl{6});
            state.ccne.phi    = y(sl{7});
            state.ccpe.phi    = y(sl{8});
            state.ccpe.E      = y(sl{9});
            
            % variables for time derivatives
            elyte_Li_cdot = yp(sl{1});
            ne_Li_csdot   = yp(sl{3});
            pe_Li_csdot   = yp(sl{5});

            elyte = model.elyte;
            ne    = model.ne;
            pe    = model.pe;
            ccne  = model.ccne;
            ccpe  = model.ccpe;

            ne_am = ne.am;
            pe_am = pe.am;

            
            elyte_cLi = state.elyte.cs{1};
            elyte_phi  = state.elyte.phi;
            ne_Li      = state.ne.am.Li;
            ne_phi     = state.ne.am.phi;
            pe_Li      = state.pe.am.Li;
            pe_phi     = state.pe.am.phi;
            ccne_phi   = state.ccne.phi;
            ccpe_phi   = state.ccpe.phi;
            ccpe_E     = state.ccpe.E;

            %% Cell voltage
            
            ccne_E = 0;
            U = ccpe_E - ccne_E;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% System of Equations                                      %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% dispatch T and SOC in submodels (need dispatch because the grids are different)

            state = model.dispatchValues(state);
            state = model.updatePhiElyte(state);
            
            state.ne.am = ne_am.updateQuantities(state.ne.am);
            state.pe.am = pe_am.updateQuantities(state.pe.am);
            
            state = setupBCSources(model, state);
            
            state.ne = ne.updateReactionRate(state.ne);
            state.pe = pe.updateReactionRate(state.pe);

            state = setupExchanges(model, state);
            
            state.elyte = elyte.updateQuantities(state.elyte);

            state.ne = ne.updateQuantities(state.ne);
            state.pe = pe.updateQuantities(state.pe);

            state.ccpe = ccpe.updateChargeCont(state.ccpe);
            state.ccne = ccne.updateChargeCont(state.ccne);
            
            %% Liquid electrolyte dissolved ionic species mass continuity and charge continuity

            elyte_Li_source = state.elyte.LiSource;
            elyte_Li_flux = state.elyte.LiFlux;

            elyte_Li_div = elyte.operators.Div(elyte_Li_flux)./elyte.G.cells.volumes;
            elyte_Li_cepsdot = elyte.eps.*elyte_Li_cdot;
            elyte_Li_massCont = (-elyte_Li_div + elyte_Li_source - elyte_Li_cepsdot);
            
            elyte_chargeCont = state.elyte.chargeCont;

            %% Electrode Active material mass continuity and charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%

            ne_Li_source = state.ne.LiSource;
            ne_Li_flux = state.ne.LiFlux;
            ne_Li_divDiff = ne.operators.Div(ne_Li_flux)./ne.G.cells.volumes;
            ne_Li_csepsdot = ne_am.eps.*ne_Li_csdot;
            ne_Li_massCont = (-ne_Li_divDiff + ne_Li_source - ne_Li_csepsdot);

            ne_e_chargeCont = state.ne.chargeCont;

            pe_Li_source = state.pe.LiSource;
            pe_Li_flux   = state.pe.LiFlux;
            pe_Li_csepsdot = pe_am.eps.*pe_Li_csdot;
            pe_Li_divDiff = model.pe.operators.Div(pe_Li_flux)./model.pe.G.cells.volumes;
            pe_Li_massCont = (-pe_Li_divDiff + pe_Li_source - pe_Li_csepsdot);

            pe_e_chargeCont =  state.pe.chargeCont;

            %% Collector charge continuity

            ccne_e_chargeCont = state.ccne.chargeCont;
            ccpe_e_chargeCont = state.ccpe.chargeCont;

            %% Control equation

            src = currentSource(t, fv.tUp, fv.tf, model.J);
            coupterm = model.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = ccpe_E;
            ccpe_sigmaeff = model.ccpe.sigmaeff;
            [tccpe, cells] = model.ccpe.operators.harmFaceBC(ccpe_sigmaeff, faces);
            control = src - sum(tccpe.*(bcval - ccpe_phi(cells)));

            %% Governing equations

            soe = vertcat(elyte_Li_massCont, ...
                          elyte_chargeCont , ...
                          ne_Li_massCont   , ...
                          ne_e_chargeCont  , ...
                          pe_Li_massCont   , ...
                          pe_e_chargeCont  , ...
                          ccne_e_chargeCont, ...
                          ccpe_e_chargeCont, ...
                          control);

        end

      function state = initStateAD(model,state)
          adbackend = model.AutoDiffBackend();
           [state.elyte.cs{1},...
           state.elyte.phi,...   
           state.ne.am.Li,...    
           state.ne.am.phi,...   
           state.pe.am.Li,...    
           state.pe.am.phi,...   
           state.ccne.phi,...    
           state.ccpe.phi,...    
           state.ccpe.E]=....
           adbackend.initVariablesAD(...
            state.elyte.cs{1},...
            state.elyte.phi,...   
            state.ne.am.Li,...    
            state.ne.am.phi,...   
            state.pe.am.Li,...    
            state.pe.am.phi,...   
            state.ccne.phi,...    
            state.ccpe.phi,...    
            state.ccpe.E);       
      end
      
     function p = getPrimaryVariables(model)
        p ={{'elyte','cs'},...
            {'elyte','phi'},...   
            {'ne','am','Li'},...    
            {'ne','am','phi'},...   
            {'pe','am','Li'},...    
            {'pe','am','phi'},...   
            {'ccne','phi'},...    
            {'ccpe','phi'},...    
            {'ccpe','E'}
          };
     end
      
      
      function state = setProp(model,state,names,val)
            nname=numel(names);
            if(nname==1)
                state.(names{1})=val;
            elseif(nname==2)
                state.(names{1}).(names{2})=val;
            elseif(nname==3)
                 state.(names{1}).(names{2}).(names{3}) = val;
            else
                error('not implmented')
            end
      end
      
      function var = getProp(model, state,names)
            var = state.(names{1});
            for i=2:numel(names)
                var = var.(names{i});
            end
      end
      
      function submod = getSubmodel(model, names)
            submod = model.(names{1});
            for i=2:numel(names)
                submod = submod.(names{i});
            end
      end
      function validforces = getValidDrivingForces(model)
         validforces=struct('src',[],'stopFunction',[]); 
      end
      
      function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
          %state_old=state;
          [problem, state] = getEquationsGen(model, state0, state,dt, drivingForces, varargin{:});
          %[pp, ss] = getEquationsGen(model, state0, state_old,dt, drivingForces, varargin{:});
          %for i=1:numel(pp.equations)
          %    assert(all(pp.equations{i}.val == problem.equations{i}.val))
          %end
              
         
      end
      
      function [problem, state] = getEquationsExp(model, state0, state,dt, drivingForces, varargin)
            time = state0.time+dt;
            % Mapping of variables
            %maybe we should have ottion to only calculate residual
            state=model.initStateAD(state);
            
            nc = model.G.cells.num;

            % setup temperature and SOC here
            %% for now this is kept constant?
            state.T =  model.T*ones(nc, 1);
            state.SOC =  model.SOC*ones(nc, 1);
            
            % variables for time derivatives
            cdotLi=struct();
            cdotLi.elyte = (state.elyte.cs{1} - state0.elyte.cs{1})/dt;
            cdotLi.ne    = (state.ne.am.Li - state0.ne.am.Li)/dt;
            cdotLi.pe   = (state.pe.am.Li - state0.pe.am.Li)/dt;

            %% Cell voltage
            
            %ccne_E = 0;
            %U = ccpe_E - ccne_E;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% System of Equations                                      %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% dispatch T and SOC in submodels (need dispatch because the grids are different)

            
            state = model.dispatchValues(state);
            state = model.updatePhiElyte(state);
            
            state.ne.am = model.ne.am.updateQuantities(state.ne.am);
            state.pe.am = model.pe.am.updateQuantities(state.pe.am);
            
            state = setupBCSources(model, state);

            state.ne = model.ne.updateReactionRate(state.ne);
            state.pe = model.pe.updateReactionRate(state.pe);
            
            state = setupExchanges(model, state);
            
            state.elyte = model.elyte.updateQuantities(state.elyte);

            state.ne = model.ne.updateQuantities(state.ne);
            state.pe = model.pe.updateQuantities(state.pe);

            state.ccpe = model.ccpe.updateChargeCont(state.ccpe);
            state.ccne = model.ccne.updateChargeCont(state.ccne);
            
            
            %% Liquid electrolyte dissolved ionic species mass continuity and charge continuity

            elyte_Li_source = state.elyte.LiSource;
            elyte_Li_flux = state.elyte.LiFlux;
            elyte_Li_div = model.elyte.operators.Div(elyte_Li_flux)./model.elyte.G.cells.volumes;
            elyte_Li_cepsdot = model.elyte.eps.*cdotLi.elyte;
            elyte_Li_massCont = (-elyte_Li_div + elyte_Li_source - elyte_Li_cepsdot);
            
            elyte_chargeCont = state.elyte.chargeCont;

            %% Electrode Active material mass continuity and charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%

            ne_Li_source = state.ne.LiSource;
            ne_Li_flux = state.ne.LiFlux;
            ne_Li_divDiff = model.ne.operators.Div(ne_Li_flux)./model.ne.G.cells.volumes;
            ne_Li_csepsdot = model.ne.am.eps.*cdotLi.ne;
            ne_Li_massCont = (-ne_Li_divDiff + ne_Li_source - ne_Li_csepsdot);

            ne_e_chargeCont = state.ne.chargeCont;

            pe_Li_source = state.pe.LiSource;
            pe_Li_flux   = state.pe.LiFlux;
            pe_Li_csepsdot = model.pe.am.eps.*cdotLi.pe;
            pe_Li_divDiff = model.pe.operators.Div(pe_Li_flux)./model.pe.G.cells.volumes;
            pe_Li_massCont = (-pe_Li_divDiff + pe_Li_source - pe_Li_csepsdot);

            pe_e_chargeCont =  state.pe.chargeCont;

            %% Collector charge continuity

            ccne_e_chargeCont = state.ccne.chargeCont;
            ccpe_e_chargeCont = state.ccpe.chargeCont;

            %% Control equation

            %src = currentSource(t, fv.tUp, fv.tf, model.J);
            src = drivingForces.src(time);%%(t, fv.tUp, fv.tf, model.J);
            coupterm = model.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = state.ccpe.E;
            ccpe_sigmaeff = model.ccpe.sigmaeff;
            [tccpe, cells] = model.ccpe.operators.harmFaceBC(ccpe_sigmaeff, faces);
            control = src - sum(tccpe.*(bcval - state.ccpe.phi(cells)));

            %% Governing equations

            eqs = {elyte_Li_massCont, ...
                   elyte_chargeCont , ...
                   ne_Li_massCont   , ...
                   ne_e_chargeCont  , ...
                   pe_Li_massCont   , ...
                   pe_e_chargeCont  , ...
                   ccne_e_chargeCont, ...
                   ccpe_e_chargeCont, ...
                   control};
               
          types={'cell','cell','cell','cell',...
              'cell','cell','cell','cell','cell'};
          names = {'elyte_Li_massCont', ...
                          'elyte_chargeCont' , ...
                          'ne_Li_massCont'   , ...
                          'ne_e_chargeCont'  , ...
                          'pe_Li_massCont'   , ...
                          'pe_e_chargeCont'  , ...
                          'ccne_e_chargeCont', ...
                          'ccpe_e_chargeCont', ...
                          'control'};
          primaryVars = model.getPrimaryVariables();            
          problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
          %state.cdotLi=cdotLi;

      end
      function [problem, state] = getEquationsGen(model, state0, state,dt, drivingForces, varargin)
              time = state0.time+dt;
              state=model.initStateAD(state);
              
              nc = model.G.cells.num;
              
              % setup temperature and SOC here
              %% for now this is kept constant?
              state.T =  model.T*ones(nc, 1);
              state.SOC =  model.SOC*ones(nc, 1);
              
              % variables for time derivatives
              cdotLi=struct();
              cdotLi.elyte = (state.elyte.cs{1} - state0.elyte.cs{1})/dt;
              cdotLi.ne    = (state.ne.am.Li - state0.ne.am.Li)/dt;
              cdotLi.pe   = (state.pe.am.Li - state0.pe.am.Li)/dt;
              
              state = model.dispatchValues(state);
              state = model.updatePhiElyte(state);
              %% first update level 2
              names={{'pe','am'},{'pe','am'}};
              for i=1:numel(names)
                  submodel = model.getSubmodel(names{i});
                  val = submodel.updateQuantities(model.getProp(state,names{i}));
                  state = model.setProp(state,names{i},val);
              end
              
              
              state = setupBCSources(model, state);
              
              names={{'ne'},{'pe'}};
              for i=1:numel(names)
                  submodel=model.getSubmodel(names{i});
                  val = submodel.updateReactionRate(model.getProp(state,names{i}));
                  state = model.setProp(state,names{i},val);
              end
              state = setupExchanges(model, state);
              
              %%update level 1
              names={{'elyte'},{'ne'},{'pe'}};
              for i=1:numel(names)
                  submodel=model.getSubmodel(names{i});
                  val = submodel.updateQuantities(model.getProp(state,names{i}));
                  state = model.setProp(state,names{i},val);
              end
              
              names={{'ccpe'},{'ccne'}};
              for i=1:numel(names)
                  submodel=model.getSubmodel(names{i});
                  val = submodel.updateChargeCont(model.getProp(state,names{i}));
                  state = model.setProp(state,names{i},val);
              end
              
              %% set equations
              names={'elyte','ne','pe'};
              eqs={};
              for i=1:numel(names)
                  submodel=model.getSubmodel({names{i}});
                  %% probably only be done on the submodel
                  source = model.getProp(state,{names{i},'LiSource'});
                  flux = model.getProp(state,{names{i},'LiFlux'});
                  %% could use submodel
                  div =  submodel.operators.Div(flux)./submodel.G.cells.volumes;
                  %% HAC
                  if(strcmp(names{i},'elyte'))
                      cepsdot = submodel.eps.*cdotLi.(names{i});
                  else
                      cepsdot = submodel.am.eps.*cdotLi.(names{i});
                  end
                  %% Li conservation
                  eqs{end+1} = -div + source - cepsdot;
                  % charge continutity
                  %% should probably be done on the sub model
                  eqs{end+1} = model.getProp(state,{names{i},'chargeCont'});
              end
              names={'ccne','ccpe'};
              for i=1:numel(names)
                  eqs{end+1} = model.getProps(state,{names{i},'chargeCont'});
              end
              
              %src = currentSource(t, fv.tUp, fv.tf, model.J);
              src = drivingForces.src(time);%%(t, fv.tUp, fv.tf, model.J);
              coupterm = model.getCoupTerm('bc-ccpe');
              faces = coupterm.couplingfaces;
              bcval = state.ccpe.E;
              ccpe_sigmaeff = model.ccpe.sigmaeff;
              [tccpe, cells] = model.ccpe.operators.harmFaceBC(ccpe_sigmaeff, faces);
              control = src - sum(tccpe.*(bcval - state.ccpe.phi(cells)));
              
              %% Governing equations
              
              eqs{end+1} = control;
          
              
              types={'cell','cell','cell','cell',...
                  'cell','cell','cell','cell','cell'};
              names = {'elyte_Li_massCont', ...
                  'elyte_chargeCont' , ...
                  'ne_Li_massCont'   , ...
                  'ne_e_chargeCont'  , ...
                  'pe_Li_massCont'   , ...
                  'pe_e_chargeCont'  , ...
                  'ccne_e_chargeCont', ...
                  'ccpe_e_chargeCont', ...
                  'control'};
              primaryVars = model.getPrimaryVariables();
              problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
              %state.cdotLi=cdotLi;
          end
          
          
          function [state, report] = updateState(model,state, problem, dx, drivingForces)
              p = model.getPrimaryVariables();
              for i=2:numel(dx)
                  val = model.getProps(state,p{i});
                  val = val + dx{i};
                  state = model.setProp(state,p{i},val);
              end
              %% not sure how to handle cells
              state.elyte.cs{1} =  state.elyte.cs{1} + dx{1};
              report = [];
          end
          
          
          
          function state = reduceState(model, state, removeContainers)
              % Reduce state to doubles, and optionally remove the property
              % containers to reduce storage space only p = model.getPrimaryVariables();
              % should be need
              state = value(state);
              
          end
          
          function [convergence, values, names] = checkConvergence(model, problem, varargin)
              
              [values, tolerances, names] = getConvergenceValues(model, problem, varargin{:});
              convergence = values < tolerances;
              fprintf('Iteration %i ',problem.iterationNo);
              fprintf(' residual ');
              fprintf(' %d ',values);
              fprintf('\n');
              disp(values)
              %disp(log(abs(values)))
          end
          
          function coupterm = getCoupTerm(model, coupname)
              coupnames = model.couplingnames;
              
              [isok, ind] = ismember(coupname, coupnames);
              assert(isok, 'name of coupling term is not recognized.');
              
              coupterm = model.couplingTerms{ind};
              
          end
          
          function coupTerm = setupNeElyteCoupTerm(model)
              
              ne = model.ne;
              elyte = model.elyte;
              
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
              
              pe = model.pe;
              elyte = model.elyte;
              
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

            ne = model.ne;
            ccne = model.ccne;

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

            pe = model.pe;
            ccpe = model.ccpe;

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

            ccne = model.ccne;
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

            ccpe = model.ccpe;
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
