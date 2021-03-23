classdef BatteryModel15i < BatteryModelSimple

    properties
        fv
    end
    
    methods

        function model = BatteryModel15i(params)
            model = model@BatteryModelSimple(params)    
        end
        
        function [t, y] = runSimulation(model)

            % Set initial conditions
            initstate = model.setupInitialState();

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

        function model = setupFV(model, state)
            model.fv = fv2d(model, state);
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

            % Here we assume E at ccne is equal to zero
            U = y(model.fv.slots{end});

            value = U - model.Ucut;
            isterminal = 1;
            direction = 0;
            
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
            elyte_phi = state.elyte.phi;
            ne_Li     = state.ne.am.Li;
            ne_phi    = state.ne.am.phi;
            pe_Li     = state.pe.am.Li;
            pe_phi    = state.pe.am.phi;
            ccne_phi  = state.ccne.phi;
            ccpe_phi  = state.ccpe.phi;
            ccpe_E    = state.ccpe.E;

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
            
            elyte_chargeCons = state.elyte.chargeCons;

            %% Electrode Active material mass continuity and charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%

            ne_Li_source = state.ne.LiSource;
            ne_Li_flux = state.ne.LiFlux;
            ne_Li_divDiff = ne.operators.Div(ne_Li_flux)./ne.G.cells.volumes;
            ne_Li_csepsdot = ne_am.eps.*ne_Li_csdot;
            ne_Li_massCont = (-ne_Li_divDiff + ne_Li_source - ne_Li_csepsdot);

            ne_e_chargeCons = state.ne.chargeCons;

            pe_Li_source = state.pe.LiSource;
            pe_Li_flux   = state.pe.LiFlux;
            pe_Li_csepsdot = pe_am.eps.*pe_Li_csdot;
            pe_Li_divDiff = model.pe.operators.Div(pe_Li_flux)./model.pe.G.cells.volumes;
            pe_Li_massCont = (-pe_Li_divDiff + pe_Li_source - pe_Li_csepsdot);

            pe_e_chargeCons =  state.pe.chargeCons;

            %% Collector charge continuity

            ccne_e_chargeCons = state.ccne.chargeCons;
            ccpe_e_chargeCons = state.ccpe.chargeCons;

            %% Control equation

            src = currentSource(t, fv.tUp, fv.tf, model.J);
            coupterm = model.getCoupTerm('bc-ccpe');
            faces = coupterm.couplingfaces;
            bcval = ccpe_E;
            ccpe_sigmaeff = model.ccpe.sigmaeff;
            [tccpe, cells] = model.ccpe.operators.harmFaceBC(ccpe_sigmaeff, faces);
            control = (src - sum(tccpe.*(bcval - ccpe_phi(cells))))/model.con.F;

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

    end

end
