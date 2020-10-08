classdef lithiumIonAD < handle
    %LITHIUMION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Identification properties
        name        % Cell name
        serialno    % Cell serial number
        form        % Cell form factor
        
        % Data structure properties
        data        % Experimental data structure
        
        % Physical constants
        con = physicalConstants();        
        
        % Domain properties
        ne          % Negative electrode object
        pe          % Positive electrode object
        sep         % Separator object
        elyte       % Electrolyte object
        
        % Model properties
        ...ecm
        ...spm
        ...p2d
        
        % Finite volume properties
        fv          % Finite volume object
        
        % State properties
        T           % Temperature,              [K]
        A           % Cell area,                [m2]
        OCV         % Cell open circuit voltage,[V]
        U           % Cell voltage,             [V]
        Ucut        % Cutoff voltage,           [V]
        I           % Cell current,             [A]
        J           % Cell current density,     [A m^-2]
        soe         % System of equations
        
        % Boundary condition properties
        chargeCont  % Galvanostatic operation BC
        
        % Preprocessing properties
        sim
        check       % Check the simulation setup, |'on'|'off'|
        
        % Postprocessing properties
        display     % Display simulation output, |'final'|'monitor'|'off'|
        style       % Style for the display

        AutoDiffBackend
    end
    
    methods
        function obj = lithiumIonAD(varargin)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin == 0
                SOC = 0.5;
                T   = 298.15;
            else
                SOC = varargin(1);
                T   = varargin(2);
            end
            
            %% Set default simulation parameters
            obj.sim{1}          = 'dynamic';        % Definition of simulation type
            obj.sim{2}          = 'isothermal';     % Definition of simulation type
            obj.check           = 'off';            % Perform user checks during preprocessing
            obj.display         = 'final';          % Regularity of output display
            obj.style.name      = 'nightfury';      % Set style of the output
            obj.style.aspect    = 1/sqrt(2);        % Set plot aspect ratio
            
            %% Define battery components
            obj.ne      = graphiteElectrode(SOC, T);
            obj.pe      = nmc111Electrode(SOC, T);
            obj.sep     = celgard2500();
            obj.elyte   = orgLiPF6(1000, T);
            
            obj.OCV     = obj.pe.am.OCP - obj.ne.am.OCP;
            obj.U       = obj.OCV;
            obj.T       = T;
            
            obj.I       = 0;
            obj.A       = 100 ./ 10000; 
            obj.J       = 1;
            obj.Ucut    = 2;
            
            obj.AutoDiffBackend = AutoDiffBackend();
        end
        
        function spm(obj)
            %SPM Creates a single particle model of the Li-ion battery
            %   Detailed explanation goes here
            
            
        end

        function [dfdy, dfdyp] = odederfun(obj, t, y, yp) 
            res = obj.odefun(t, y, yp, 'useAD', true);
            dfdy = res.jac{1};
            dfdyp = res.jac{2};
        end
    
         
        function [t, y] = p2d(obj)
           
            % Generate the FV mesh
            names           = {'NE';'SEP';'PE'};
            sizes           = [ 1e-6; 1e-6; 1e-6 ];
            lengths         = [ obj.ne.t; ...
                                obj.sep.t; ...
                                obj.pe.t ];
                            
            obj.fv = fv1d(names, sizes, lengths);
            
            obj.ne.dombin   = obj.fv.dombin{1};
            obj.sep.dombin  = obj.fv.dombin{2};
            obj.pe.dombin   = obj.fv.dombin{3};
            
            obj.ne.N        = obj.fv.Nvec(1);
            obj.sep.N       = obj.fv.Nvec(2);
            obj.pe.N        = obj.fv.Nvec(3);
            obj.elyte.N     = obj.fv.N;
            
            obj.ne.X        = obj.fv.X(obj.ne.dombin>0);
            obj.sep.X       = obj.fv.X(obj.sep.dombin>0);
            obj.pe.X        = obj.fv.X(obj.pe.dombin>0);
            obj.elyte.X     = obj.fv.X;
            
            obj.ne.Xb       = obj.fv.Xb(1:obj.ne.N+1);
            obj.sep.Xb      = obj.fv.Xb(obj.ne.N+1:obj.ne.N + obj.sep.N + 1);
            obj.pe.Xb       = obj.fv.Xb(obj.ne.N + obj.sep.N + 1:end);
            obj.elyte.Xb    = obj.fv.Xb;
            
            obj.icp2d();
            ff={'ne','sep','pe','elyte'};
            for i=1:numel(ff)
                dd=ff{i};
                obj.(dd).Grid = cartGrid([obj.(dd).N],[max(obj.(dd).Xb)-min(obj.(dd).Xb)]);
                obj.(dd).Grid.nodes.coords(:,1)=obj.(dd).Grid.nodes.coords(:,1)+min(obj.(dd).Xb);
                obj.(dd).Grid = computeGeometry(obj.(dd).Grid);
                rock = struct('perm',ones(obj.(dd).Grid.cells.num,1),'poro',ones(obj.(dd).Grid.cells.num,1));
                obj.(dd).operators = setupOperatorsTPFA(obj.(dd).Grid,rock);
                obj.(dd).operators.allDiv = getAllDiv(obj.(dd).Grid)
            end
            % Set initial conditions
            obj.icp2d();
            %% Is this constant?obj.elyte.sp.Li.Deff
            [~,obj.elyte.sp.Li.Trans] = divAgradCreate( obj.elyte.sp.Li.c, ...
                                                    -obj.elyte.sp.Li.Deff, ...
                                                    obj.fv.X, ...
                                                    obj.fv.Xb);  
            
            [~,obj.ne.am.Li.Trans] = divAgradCreate( obj.ne.am.Li.cs, ...
                                             -obj.ne.am.Li.Deff, ...
                                             obj.ne.X, ...
                                             obj.ne.Xb);
            
            [~,obj.pe.am.Li.Trans] = divAgradCreate( obj.pe.am.Li.cs, ...
                                             -obj.pe.am.Li.Deff, ...
                                             obj.pe.X, ...
                                             obj.pe.Xb);
            % Perform manual mesh chec
            if strcmpi(obj.check, 'on')==1
                obj.figureWindow()
                obj.plotMesh()
                ui = input('Is the mesh ok? (yes = 1, no = 0)    : ');
                if ui == 0
                    error('Execution stopped by user: mesh is not ok')
                else
                    close all
                end
            end
            
           
        
            %% Solution space (time or current) discretization
                
            % Time discretization
            obj.fv.ti = 0;
            obj.fv.tf = 3600*24;
            obj.fv.tf = 30;
            obj.fv.dt = 10;
            obj.fv.tUp = 0.1;
            obj.fv.tSpan = obj.fv.ti:obj.fv.dt:obj.fv.tf;

            % Pre-process
            obj.dynamicPreprocess();

            % Solve the system of equations
            endFun = @(t,y,yp)obj.cutOff(t,y,yp);
            
            fun = @(t,y,yp) obj.odefun(t, y, yp);
            derfun = @(t,y,yp) obj.odederfun(t, y, yp);           

            options = odeset('RelTol'  , 1e-4  , ...
                             'AbsTol'  , 1e-6  , ... 
                             'Stats'   , 'on'  , ... 
                             'Events'  , endFun, ...
                             'Jacobian', derfun);
            
            [t, y] = ode15i(fun, obj.fv.tSpan', obj.fv.y0, obj.fv.yp0, options);
                
            doplot = true;
            if doplot
                resultname = 'test';
                path = 'C:\Users\simonc\Documents\01_Projects\Electrochemical Modelling Toolbox\Batteries\';
                day = date();
                hrmin = clock();
                hr = num2str(hrmin(4));
                minhr = num2str(hrmin(5));
                save([path, resultname, '-', day,'-', hr,'h', minhr,'m','.mat']);
                
                if strcmpi(obj.display, 'final')
                    obj.plotSummary(t,y);
                end
            end
            
        end
        
        function [value,isterminal,direction] = cutOff(obj,t,y,yp)
            %UNTITLED3 Summary of this function goes here
            %   Detailed explanation goes here

            %value = sum(Zn_s.Asp.*Netz.deltaX) - 0.2*sum(Zn_s.Asp_0.*Netz.deltaX);

            value = (obj.U) - obj.Ucut;
            isterminal = 1;
            direction = 0;

        end
        
        function icp2d(obj)
            
            %% Negative electrode properties
            % Temperature
            obj.ne.am.T             = obj.ne.am.T .* ones(obj.ne.N, 1);
            
            % Volume fractions
                % Solid volume fractions
            obj.ne.am.eps           = obj.ne.am.eps  .* ones(obj.ne.N, 1);
            obj.ne.bin.eps          = obj.ne.bin.eps .* ones(obj.ne.N, 1);
            obj.ne.sei.eps          = obj.ne.sei.eps .* ones(obj.ne.N, 1);
            obj.ne.eps              = obj.ne.am.eps ...
                                    + obj.ne.bin.eps ...
                                    + obj.ne.sei.eps;
                % Void volume fraction
            obj.ne.void             = 1 - obj.ne.eps;
                            
            % SEI thickness
            obj.ne.sei.t            = obj.ne.sei.t .* ones(obj.ne.N, 1);
            
            % Lithiation
            obj.ne.am.Li.cs         = obj.ne.am.Li.cs .* ones(obj.ne.N, 1);
            obj.ne.am.Li.cseps      = obj.ne.am.Li.cs .* obj.ne.am.eps;
            obj.ne.am.theta         = obj.ne.am.theta .* ones(obj.ne.N, 1);
            obj.ne.am.SOC           = obj.ne.am.SOC   .* ones(obj.ne.N, 1);
            
            % Open-circuit potential
            obj.ne.am.OCP           = obj.ne.am.OCP .* ones(obj.ne.N, 1);
            obj.ne.am.phi           = obj.ne.am.OCP;
            
            % Kinetic Coefficients
            obj.ne.am.k             = obj.ne.am.k .* ones(obj.ne.N, 1);
            
            % Transport Coefficients
            obj.ne.am.Li.D          = obj.ne.am.Li.D .* ones(obj.ne.N, 1);
            obj.ne.am.Li.Deff       = obj.ne.am.Li.D .* obj.ne.am.eps.^1.5;
            
            %% Positive electrode properties
            % Temperature
            obj.pe.am.T             = obj.pe.am.T .* ones(obj.pe.N, 1);
            
            % Volume fractions
                % Solid volume fractions
            obj.pe.am.eps           = obj.pe.am.eps .* ones(obj.pe.N, 1);
            obj.pe.bin.eps          = obj.pe.bin.eps .* ones(obj.pe.N, 1);
            obj.pe.eps              = obj.pe.am.eps ...
                                    + obj.pe.bin.eps;
                % Void volume fraction
            obj.pe.void             = 1 - obj.pe.eps;
            
            % Lithiation
            obj.pe.am.Li.cs         = obj.pe.am.Li.cs .* ones(obj.pe.N, 1);
            obj.pe.am.Li.cseps      = obj.pe.am.Li.cs .* obj.pe.am.eps;
            obj.pe.am.theta         = obj.pe.am.theta .* ones(obj.pe.N, 1);
            obj.pe.am.SOC           = obj.pe.am.SOC .* ones(obj.pe.N, 1);
            
            % Open-circuit potential
            obj.pe.am.OCP           = obj.pe.am.OCP .* ones(obj.pe.N, 1);
            obj.pe.am.phi           = obj.pe.am.OCP;
            
            % Kinetic Coefficients
            obj.pe.am.k             = obj.pe.am.k .* ones(obj.pe.N, 1);
            
            % Transport Coefficients
            obj.pe.am.Li.D          = obj.pe.am.Li.D .* ones(obj.pe.N, 1);
            obj.pe.am.Li.Deff       = obj.pe.am.Li.D .* obj.pe.am.eps.^1.5;
            
            %% Separator properties
            obj.sep.eps             = obj.sep.eps .* ones(obj.sep.N, 1);
            obj.sep.void            = 1 - obj.sep.eps;
            
            %% Electrolyte properties
            obj.elyte.eps           = [ obj.ne.void; ...
                                        obj.sep.void; ...
                                        obj.pe.void ];
                            
            obj.elyte.phi           = zeros(obj.elyte.N, 1);
            obj.elyte.sp.Li.c       = obj.elyte.sp.Li.c .* ones(obj.elyte.N, 1);
            obj.elyte.sp.Li.ceps    = obj.elyte.sp.Li.c .* obj.elyte.eps;
            obj.elyte.T             = obj.elyte.T .* ones(obj.elyte.N, 1);
            obj.elyte.sp.Li.D       = obj.elyte.sp.Li.D .* ones(obj.elyte.N, 1);
            obj.elyte.sp.Li.Deff    = obj.elyte.sp.Li.D .* obj.elyte.eps .^1.5;
            obj.elyte.kappa         = obj.elyte.kappa .* ones(obj.elyte.N, 1);
            obj.elyte.kappaeff      = obj.elyte.kappa .* obj.elyte.eps .^1.5;
            
        end
        
        function dynamicPreprocess(obj)
            
            %% Initialize the state vector
            obj.fv.y0   = [ obj.elyte.sp.Li.ceps;
                            obj.elyte.phi;
                            obj.ne.am.Li.cseps;
                            obj.ne.am.phi;
                            obj.pe.am.Li.cseps;
                            obj.pe.am.phi;
                            obj.pe.E];
            
            %% Inititalize the state time derivative vector
            obj.fv.yp0  = zeros(length(obj.fv.y0),1);
            
            %% Store state slots
            obj.fv.s1   = 1:obj.elyte.N; 
            obj.fv.s2   = (obj.fv.s1(end)+1):(obj.fv.s1(end)+obj.elyte.N);
            obj.fv.s3   = (obj.fv.s2(end)+1):(obj.fv.s2(end)+obj.ne.N);
            obj.fv.s4   = (obj.fv.s3(end)+1):(obj.fv.s3(end)+obj.ne.N);
            obj.fv.s5   = (obj.fv.s4(end)+1):(obj.fv.s4(end)+obj.pe.N);
            obj.fv.s6   = (obj.fv.s5(end)+1):(obj.fv.s5(end)+obj.pe.N);
            obj.fv.s7   = length(obj.fv.y0);
            
        end
        
        function dynamicReadState(obj, y, yp, varargin)
            
            opt = struct('useAD', false);
            opt = merge_options(opt, varargin{:});
            useAD = opt.useAD;
            
            if useAD
                adbackend = obj.AutoDiffBackend();
                [y, yp] = adbackend.initVariablesAD(y, yp);
            end
            
            obj.elyte.sp.Li.ceps    = y(obj.fv.s1);
            obj.elyte.phi           = y(obj.fv.s2);
            obj.ne.am.Li.cseps      = y(obj.fv.s3);
            obj.ne.am.phi           = y(obj.fv.s4);
            obj.pe.am.Li.cseps      = y(obj.fv.s5);
            obj.pe.am.phi           = y(obj.fv.s6);
            obj.pe.E                = y(end);
            
            obj.elyte.sp.Li.cepsdot = yp(obj.fv.s1);
            obj.ne.am.Li.csepsdot   = yp(obj.fv.s3);
            obj.pe.am.Li.csepsdot   = yp(obj.fv.s5);
            
            obj.elyte.sp.Li.c  = obj.elyte.sp.Li.ceps ./ obj.elyte.eps;
            obj.elyte.sp.PF6.c = obj.elyte.sp.Li.c;
            
            obj.ne.am.Li.cs = obj.ne.am.Li.cseps ./ obj.ne.am.eps;
            obj.pe.am.Li.cs = obj.pe.am.Li.cseps ./ obj.pe.am.eps;
            
            %% Update electrolyte physicochemical and transport properties
            obj.elyte.update()
            obj.ne.am.update()
            obj.pe.am.update()
            
            obj.elyte.kappaeff = obj.elyte.kappa .* ...
                obj.elyte.eps .^1.5;
            obj.elyte.sp.Li.Deff = obj.elyte.sp.Li.D .* ...
                obj.elyte.eps .^1.5;
            
            obj.ne.am.Li.Deff  = obj.ne.am.Li.D .* obj.ne.am.eps.^1.5;
            obj.ne.sigmaeff = obj.ne.am.sigma .* obj.ne.am.eps.^1.5;
            obj.pe.am.Li.Deff  = obj.pe.am.Li.D .* obj.pe.am.eps.^1.5;
            obj.pe.sigmaeff = obj.pe.am.sigma .* obj.pe.am.eps.^1.5;
            
            
            %% Ionic current density   
            % Ionic current density in the liquid
            %   Ionic current density due to the chemical potential gradient
            N = obj.elyte.N;
            ncomp = obj.elyte.ncomp;
            
            jchems = cell(ncomp, 1);
            for i = 1 : ncomp
                jchems{i} = zeros(N + 1, 1);         
                jchems{i} = harm(obj.elyte.kappaeff .* obj.elyte.ion.tvec{i} .* obj.elyte.ion.dmudc{i} ./ ...
                                 (obj.elyte.ion.zvec{i}.*obj.con.F), obj.fv.X, obj.fv.Xb) .* grad(obj.elyte.ion.cvec{i}, ...
                                                                  obj.fv.X);
            end
            obj.elyte.jchem = jchems{1};
            for i = 2 : ncomp
                obj.elyte.jchem = obj.elyte.jchem + jchems{i};
            
            end
            %   Ionic current density due to the electrochemical potential gradient
            obj.elyte.j = harm(obj.elyte.kappaeff, obj.fv.X, obj.fv.Xb).*(-1).*grad(obj.elyte.phi, obj.fv.X) ...
                - obj.elyte.jchem;   
            
            %% Electric current density
            % Active material NE
            obj.ne.j =  harm(obj.ne.sigmaeff, obj.ne.X, obj.ne.Xb) ...
                        .* (-1) .* grad(obj.ne.am.phi, ...
                                        obj.ne.X, ...
                                        'left', ...
                                        'dirlichet', ...
                                        0);
                                    
            % Active material PE                        
            obj.pe.j = harm(obj.pe.sigmaeff, obj.pe.X, obj.pe.Xb) ...
                        .* (-1) .* grad(obj.pe.am.phi, ...
                                        obj.pe.X, ...
                                        'right',...
                                        'dirlichet', ...
                                        obj.pe.E);
                                    
            %% Cell voltage
            obj.U = obj.pe.E - obj.ne.E;
        end
        
        function dynamicBuildSOE(obj, t, varargin)
            
            opt = struct('useAD', false);
            opt = merge_options(opt, varargin{:});
            useAD = opt.useAD;
            
            if useAD
                adsample = getSampleAD(obj.elyte.sp.Li.ceps, ...    
                                       obj.elyte.phi, ...           
                                       obj.ne.am.Li.cseps, ...      
                                       obj.ne.am.phi, ...           
                                       obj.pe.am.Li.cseps, ...      
                                       obj.pe.am.phi, ...           
                                       obj.pe.E, ...                
                                       obj.elyte.sp.Li.cepsdot, ... 
                                       obj.ne.am.Li.csepsdot, ...   
                                       obj.pe.am.Li.csepsdot);
                adbackend = obj.AutoDiffBackend;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Source terms for continuity equations                    %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Inititalize vectors
            obj.elyte.sp.Li.source  = zeros(obj.fv.N, 1);
            obj.ne.am.Li.source     = zeros(obj.ne.N, 1);
            obj.pe.am.Li.source     = zeros(obj.pe.N, 1);
            obj.ne.am.e.source      = zeros(obj.ne.N, 1);
            obj.pe.am.e.source      = zeros(obj.pe.N, 1);
            
            if useAD
                obj.elyte.sp.Li.source = adbackend.convertToAD(obj.elyte.sp.Li.source, adsample);
                obj.ne.am.Li.source    = adbackend.convertToAD(obj.ne.am.Li.source, adsample);
                obj.pe.am.Li.source    = adbackend.convertToAD(obj.pe.am.Li.source, adsample);
                obj.ne.am.e.source     = adbackend.convertToAD(obj.ne.am.e.source, adsample);
                obj.pe.am.e.source     = adbackend.convertToAD(obj.pe.am.e.source, adsample);
            end
            
            % Calculate reaction rates           
            obj.ne.reactBV(obj.elyte.phi(1:obj.ne.N));
            obj.pe.reactBV(obj.elyte.phi(obj.ne.N + obj.sep.N + 1:end));
            
            % Set source terms            
            %%%%% Li+ Sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Electrolyte NE Li+ source
            obj.elyte.sp.Li.source(1:obj.ne.N) = +1 .* obj.ne.R;

            % Electrolyte PE Li+ source
            obj.elyte.sp.Li.source(obj.ne.N+obj.sep.N+1:end) = -1 .* obj.pe.R;
            
            %%%%% Li0 Sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Active Material NE Li0 source
            obj.ne.am.Li.source = -1 .* obj.ne.R;

            % Active Material PE Li0 source
            obj.pe.am.Li.source = +1 .* obj.pe.R;
            
            %%%%% e- Sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Active Material NE Li0 source
            obj.ne.am.e.source = +1 .* obj.ne.R;

            % Active Material PE Li0 source
            obj.pe.am.e.source = -1 .* obj.pe.R;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Diffusion Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of diffusion mass flux
            %   Electrolyte Li+ Diffusion
            if(false)
            obj.elyte.sp.Li.divDiff = divAgradCreate( obj.elyte.sp.Li.c, ...
                                                -obj.elyte.sp.Li.Deff, ...
                                                obj.fv.X, ...
                                                obj.fv.Xb);  
            obj.ne.am.Li.divDiff = divAgradCreate( obj.ne.am.Li.cs, ...
                                             -obj.ne.am.Li.Deff, ...
                                             obj.ne.X, ...
                                             obj.ne.Xb);
            
            obj.pe.am.Li.divDiff = divAgradCreate( obj.pe.am.Li.cs, ...
                                             -obj.pe.am.Li.Deff, ...
                                             obj.pe.X, ...
                                             obj.pe.Xb);
            else
               assert( isa(obj.pe.am.Li.Deff,'double') && isa(obj.ne.am.Li.Deff,'double') && isa(obj.elyte.sp.Li.Deff,'double') )
            x = obj.elyte.sp.Li.c;
            flux = obj.elyte.sp.Li.Trans.*obj.elyte.operators.Grad(x);
            obj.elyte.sp.Li.divDiff=  obj.elyte.operators.Div(flux)./obj.elyte.Grid.cells.volumes;
            %%
            x= obj.ne.am.Li.cs;
            flux = obj.ne.am.Li.Trans.*obj.ne.operators.Grad(x);
            obj.ne.am.Li.divDiff=  obj.ne.operators.Div(flux)./obj.ne.Grid.cells.volumes;
            
            %%
            x = obj.pe.am.Li.cs;
            flux = obj.pe.am.Li.Trans.*obj.pe.operators.Grad(x);
            obj.pe.am.Li.divDiff=  obj.pe.operators.Div(flux)./obj.pe.Grid.cells.volumes;
            
            
                
                
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Migration Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of the migration mass flux
            %   Electrolyte Li+ Migration
            flux = obj.elyte.sp.Li.t ./ (obj.elyte.sp.Li.z .* obj.con.F) ...
                                                .* obj.elyte.j;
            %obj.elyte.sp.Li.divMig    = div( flux, obj.fv.Xb);
                                            
            obj.elyte.sp.Li.divMig   =  obj.elyte.operators.allDiv(flux);%./obj.elyte.Grid.cells.volumes;                           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% System of Equations                                      %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Liquid electrolyte dissolved ionic species mass continuity %
            %   Electrolyte Li+ Mass Continuity
            obj.elyte.sp.Li.massCont =    (   - obj.elyte.sp.Li.divDiff ...
                                                - obj.elyte.sp.Li.divMig  ...
                                                + obj.elyte.sp.Li.source  ...
                                                - obj.elyte.sp.Li.cepsdot);
           
            %% Liquid electrolyte charge continuity %%%%%%%%%%%%%%%%%%%%%%%
            
            %obj.elyte.chargeCont = -div(  obj.elyte.j, obj.fv.Xb) ./ obj.con.F + ...
            %                         obj.elyte.sp.Li.source .* obj.elyte.sp.Li.z;       
            obj.elyte.chargeCont = -obj.elyte.operators.allDiv(  obj.elyte.j) ./ obj.con.F + ...
                                     obj.elyte.sp.Li.source .* obj.elyte.sp.Li.z;
            %a = -div(  obj.elyte.j, obj.fv.Xb);
            %b = -obj.elyte.operators.allDiv(  obj.elyte.j)
            %assert(all(a==b))
             %% Active material mass continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%
             obj.ne.am.Li.massCont = (-obj.ne.am.Li.divDiff ...
                                        + obj.ne.am.Li.source ...
                                        - obj.ne.am.Li.csepsdot);
                                    
             obj.pe.am.Li.massCont = (-obj.pe.am.Li.divDiff ...
                                            + obj.pe.am.Li.source ...
                                            - obj.pe.am.Li.csepsdot);
                                        
            %% Active material charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%
            % obj.ne.am.e.chargeCont = -div(  obj.ne.j, obj.ne.Xb) ./ -obj.con.F + ...
            %                         obj.ne.am.e.source .* -1;   
            obj.ne.am.e.chargeCont = -obj.ne.operators.allDiv(  obj.ne.j) ./ -obj.con.F + ...
                                     obj.ne.am.e.source .* -1;
                                 
           %obj.pe.am.e.chargeCont = -div(  obj.pe.j, obj.pe.Xb) ./ -obj.con.F + ...
           %                          obj.pe.am.e.source .* -1;   
           obj.pe.am.e.chargeCont = -obj.pe.operators.allDiv(  obj.pe.j) ./ -obj.con.F + ...
                                     obj.pe.am.e.source .* -1;              
            %% Global charge continuity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.chargeCont = currentSource(t, obj.fv.tUp, obj.fv.tf, obj.J) ...
                            - (sum(obj.pe.R .* obj.con.F .* obj.fv.dXvec(3)));
            
            %% State vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
            obj.soe = vertcat(obj.elyte.sp.Li.massCont, ...
                              obj.elyte.chargeCont, ...
                              obj.ne.am.Li.massCont, ...
                              obj.ne.am.e.chargeCont, ...
                              obj.pe.am.Li.massCont, ...
                              obj.pe.am.e.chargeCont, ...
                              obj.chargeCont);

        end
        
        %function a = divAgrad(y)
            
        %end
        
        function res = odefun(obj, t, y, yp, varargin)
        %ODEFUN Compiles the system of differential equations
            
            opt = struct('useAD', false);
            opt = merge_options(opt, varargin{:});
            useAD = opt.useAD;
            
            dodisplay = false;
            if dodisplay
                % Display current time step
                clc
                disp('Time, s:')
                disp(t)
                disp('Current Density, A/cm2:')
                disp(-1e-4.*currentSource(t, obj.fv.tUp, obj.fv.tf, obj.J))
                disp('Cell Voltage, V:')
                disp(obj.pe.E - obj.ne.E)
            end
            
            % Read the state vector
            obj.dynamicReadState(y, yp, 'useAD', useAD);

            % Build SOE
            obj.dynamicBuildSOE(t, 'useAD', useAD);
            
            res = obj.soe;
        end     
        
  
        

        
        
        function plotMesh(obj)
            
            subplot(3,1,1), plot(obj.fv.Xb.*1e3, zeros(size(obj.fv.Xb)),'+','MarkerSize',5, 'LineWidth',2, 'Color', obj.style.palette.discrete(1,:));
            hold on
            plot(obj.fv.X.*1e3, zeros(size(obj.fv.X)),'o','MarkerSize',5, 'LineWidth',2, 'Color', obj.style.palette.discrete(7,:));
            hold off
            xlim([min(obj.fv.Xb), obj.fv.Xb(10)].*1e3);
            set(gca,'FontSize',28, ...
                'color',obj.style.background, ...
                'ColorOrder', obj.style.palette.discrete)
            ax = gca;
            ax.XColor = 'w';
            ax.YColor = 'w';
            title('Finite Volume Mesh: Left Bound', 'Color', 'w')
            
            subplot(3,1,2), plot(obj.fv.Xb.*1e3, zeros(size(obj.fv.Xb)),'+','MarkerSize',5, 'LineWidth',2, 'Color', obj.style.palette.discrete(1,:));
            hold on
            plot(obj.fv.X.*1e3, zeros(size(obj.fv.X)),'o','MarkerSize',5, 'LineWidth',2, 'Color', obj.style.palette.discrete(7,:));
            hold off
            xlim([min(obj.fv.Xb), max(obj.fv.Xb)].*1e3);
            set(gca,'FontSize',28, ...
                'color',obj.style.background, ...
                'ColorOrder', obj.style.palette.discrete)
            ax = gca;
            ax.XColor = 'w';
            ax.YColor = 'w';
            title('Finite Volume Mesh: Whole Domain', 'Color', 'w')
            
            subplot(3,1,3), plot(obj.fv.Xb.*1e3, zeros(size(obj.fv.Xb)),'+','MarkerSize',5, 'LineWidth',2, 'Color', obj.style.palette.discrete(1,:));
            hold on
            plot(obj.fv.X.*1e3, zeros(size(obj.fv.X)),'o','MarkerSize',5, 'LineWidth',2, 'Color', obj.style.palette.discrete(7,:));
            hold off
            xlim([obj.fv.Xb(end-10), max(obj.fv.Xb)].*1e3);
            set(gca,'FontSize',28, ...
                'color',obj.style.background, ...
                'ColorOrder', obj.style.palette.discrete)
            ax = gca;
            ax.XColor = 'w';
            ax.YColor = 'w';
            title('Finite Volume Mesh: Right Bound', 'Color', 'w')
            hold off
            xlabel('Location  /  mm')
            
        end
        
        function plotElyteConcentration(obj, varargin)
            
            if isempty(varargin)
                plot(obj.fv.X.*1e3, obj.elyte.sp.Li.c.*1e-3, 'LineWidth', 5)
                xlabel('Position  /  mm')
                ylabel('LiPF_6 Conc.  /  mol\cdotL^{-1}')
                eval(obj.style.name)
            end
            
        end
        
        function plotElytePotential(obj, varargin)
            
            if isempty(varargin)
                plot(obj.fv.X.*1e3, obj.elyte.phi, 'LineWidth', 5)
                xlabel('Position  /  mm')
                ylabel('Electric Potential  /  V')
                eval(obj.style.name)
            end
            
        end
        
        function plotPotentials(obj, varargin)
            
            if isempty(varargin)
                plot(obj.fv.X.*1e3, obj.elyte.phi, 'LineWidth', 5)
                hold on
                plot(obj.ne.X.*1e3, obj.ne.am.phi, 'LineWidth', 5)
                plot(obj.pe.X.*1e3, obj.pe.am.phi, 'LineWidth', 5)
                xlabel('Position  /  mm')
                ylabel('Electric Potential  /  V')
                eval(obj.style.name)
            end
            
        end
        
        function plotCurrentDensity(obj, varargin)
            
            if isempty(varargin)
                plot(obj.fv.Xb.*1e3, obj.elyte.j, 'LineWidth', 5)
                hold on
                plot(obj.ne.Xb.*1e3, obj.ne.j, 'LineWidth', 5)
                plot(obj.pe.Xb.*1e3, obj.pe.j, 'LineWidth', 5)
                xlabel('Position  /  mm')
                ylabel('Current Density  /  A m^{-2}')
                eval(obj.style.name)
            end
            
        end
        
        function plotSummary(obj,t,y,varargin)
            
            close all
            
            if ~isempty(varargin)
                teval = varargin;
                [~,~,idt]=unique(round(abs(t-teval)),'stable');
            else
                idt = length(t);
            end
            
          figure(1), plot(t./3600, y(:,end)-obj.ne.E, 'LineWidth', 5) 
           ylim([2, 4.2])
           xlabel('Time  /  h')
           ylabel('Cell Voltage  /  V')
           eval(obj.style.name)
           
           figure(2)
           subplot(2,3,1:3), plot(obj.fv.X.*1e3, 1e-3.*y(idt, obj.fv.s1)./obj.elyte.eps', 'LineWidth', 5)
                xlabel('Position  /  mm')
                ylabel('LiPF_6 Conc.  /  mol\cdotL^{-1}')
                eval([obj.style.name, '("subplot")'])
                set(gca, 'FontSize', 14)
                
            subplot(2,3,4), plot(obj.ne.X.*1e3, y(idt, obj.fv.s3)./obj.ne.am.eps' ./ obj.ne.am.Li.cmax, 'LineWidth', 5)
                xlabel('Position  /  mm')
                %ylabel('Li Conc.  /  mol\cdotL^{-1}')
                eval([obj.style.name, '("subplot")'])
                set(gca, 'FontSize', 14)
                
            subplot(2,3,6), plot(obj.pe.X.*1e3, y(idt, obj.fv.s5)./obj.pe.am.eps'./ obj.pe.am.Li.cmax, 'LineWidth', 5)
                xlabel('Position  /  mm')
                %ylabel('Li Conc.  /  mol\cdotL^{-1}')
                eval([obj.style.name, '("subplot")'])
                set(gca, 'FontSize', 14)
                
            figure(3)
           subplot(2,3,1:3), plot(obj.fv.X.*1e3, y(idt, obj.fv.s2), 'LineWidth', 5)
                xlabel('Position  /  mm')
                ylabel('Electric Potential  /  V')
                eval([obj.style.name, '("subplot")'])
                set(gca, 'FontSize', 14)
                
            subplot(2,3,4), plot(obj.ne.X.*1e3, y(idt, obj.fv.s4), 'LineWidth', 5)
                xlabel('Position  /  mm')
                %ylabel('Li Conc.  /  mol\cdotL^{-1}')
                eval([obj.style.name, '("subplot")'])
                set(gca, 'FontSize', 14)
                
            subplot(2,3,6), plot(obj.pe.X.*1e3, y(idt, obj.fv.s6), 'LineWidth', 5)
                xlabel('Position  /  mm')
                %ylabel('Li Conc.  /  mol\cdotL^{-1}')
                eval([obj.style.name, '("subplot")'])
                set(gca, 'FontSize', 14)
           
            
        end
        
        function figureWindow(obj,font)
            figure()
            set(gca,'FontSize',28, ...
                'color',obj.style.background, ...
                'ColorOrder', obj.style.palette.discrete)
            ax = gca;
            ax.XColor = 'w';
            ax.YColor = 'w';
            set(gcf,'units','centimeter',...
                'position',[5,5,obj.style.width,obj.style.height], ...
                'color',obj.style.background)
            if nargin < 2 
                set(gca, 'FontName', 'Helvetica')
            elseif srtcmpi(font, 'serif') == 1
                set(gca, 'FontName', 'Baskerville Old Face')
            else
                set(gca, 'FontName', 'Helvetica')
            end     
            hold on
        end
        
    end
    
    methods(Static)
        function map = magma(~,N)
            %MAGMA Defines the properties of the magma colormap
            %   obj.magma() calculates the magma colormap           
            
            if nargin<2 || isempty(N)
                N = 256;%size(get(gcf,'colormap'),1);
            else
                assert(isscalar(N)&&isreal(N),'First argument must be a real numeric scalar.')
            end
            
            % Define the magma colormap
            raw = [0.001462,0.000466,0.013866; 
                0.002258,0.001295,0.018331; 
                0.003279,0.002305,0.023708; 
                0.004512,0.003490,0.029965; 
                0.005950,0.004843,0.037130; 
                0.007588,0.006356,0.044973; 
                0.009426,0.008022,0.052844; 
                0.011465,0.009828,0.060750; 
                0.013708,0.011771,0.068667; 
                0.016156,0.013840,0.076603; 
                0.018815,0.016026,0.084584; 
                0.021692,0.018320,0.092610; 
                0.024792,0.020715,0.100676; 
                0.028123,0.023201,0.108787; 
                0.031696,0.025765,0.116965; 
                0.035520,0.028397,0.125209; 
                0.039608,0.031090,0.133515; 
                0.043830,0.033830,0.141886; 
                0.048062,0.036607,0.150327; 
                0.052320,0.039407,0.158841; 
                0.056615,0.042160,0.167446; 
                0.060949,0.044794,0.176129; 
                0.065330,0.047318,0.184892; 
                0.069764,0.049726,0.193735; 
                0.074257,0.052017,0.202660; 
                0.078815,0.054184,0.211667; 
                0.083446,0.056225,0.220755; 
                0.088155,0.058133,0.229922; 
                0.092949,0.059904,0.239164; 
                0.097833,0.061531,0.248477; 
                0.102815,0.063010,0.257854; 
                0.107899,0.064335,0.267289; 
                0.113094,0.065492,0.276784; 
                0.118405,0.066479,0.286321; 
                0.123833,0.067295,0.295879; 
                0.129380,0.067935,0.305443; 
                0.135053,0.068391,0.315000; 
                0.140858,0.068654,0.324538; 
                0.146785,0.068738,0.334011; 
                0.152839,0.068637,0.343404; 
                0.159018,0.068354,0.352688; 
                0.165308,0.067911,0.361816; 
                0.171713,0.067305,0.370771; 
                0.178212,0.066576,0.379497; 
                0.184801,0.065732,0.387973; 
                0.191460,0.064818,0.396152; 
                0.198177,0.063862,0.404009; 
                0.204935,0.062907,0.411514; 
                0.211718,0.061992,0.418647; 
                0.218512,0.061158,0.425392; 
                0.225302,0.060445,0.431742; 
                0.232077,0.059889,0.437695; 
                0.238826,0.059517,0.443256; 
                0.245543,0.059352,0.448436; 
                0.252220,0.059415,0.453248; 
                0.258857,0.059706,0.457710; 
                0.265447,0.060237,0.461840; 
                0.271994,0.060994,0.465660; 
                0.278493,0.061978,0.469190; 
                0.284951,0.063168,0.472451; 
                0.291366,0.064553,0.475462; 
                0.297740,0.066117,0.478243; 
                0.304081,0.067835,0.480812; 
                0.310382,0.069702,0.483186; 
                0.316654,0.071690,0.485380; 
                0.322899,0.073782,0.487408; 
                0.329114,0.075972,0.489287; 
                0.335308,0.078236,0.491024; 
                0.341482,0.080564,0.492631; 
                0.347636,0.082946,0.494121; 
                0.353773,0.085373,0.495501; 
                0.359898,0.087831,0.496778; 
                0.366012,0.090314,0.497960; 
                0.372116,0.092816,0.499053; 
                0.378211,0.095332,0.500067; 
                0.384299,0.097855,0.501002; 
                0.390384,0.100379,0.501864; 
                0.396467,0.102902,0.502658; 
                0.402548,0.105420,0.503386; 
                0.408629,0.107930,0.504052; 
                0.414709,0.110431,0.504662; 
                0.420791,0.112920,0.505215; 
                0.426877,0.115395,0.505714; 
                0.432967,0.117855,0.506160; 
                0.439062,0.120298,0.506555; 
                0.445163,0.122724,0.506901; 
                0.451271,0.125132,0.507198; 
                0.457386,0.127522,0.507448; 
                0.463508,0.129893,0.507652; 
                0.469640,0.132245,0.507809; 
                0.475780,0.134577,0.507921; 
                0.481929,0.136891,0.507989; 
                0.488088,0.139186,0.508011; 
                0.494258,0.141462,0.507988; 
                0.500438,0.143719,0.507920; 
                0.506629,0.145958,0.507806; 
                0.512831,0.148179,0.507648; 
                0.519045,0.150383,0.507443; 
                0.525270,0.152569,0.507192; 
                0.531507,0.154739,0.506895; 
                0.537755,0.156894,0.506551; 
                0.544015,0.159033,0.506159; 
                0.550287,0.161158,0.505719; 
                0.556571,0.163269,0.505230; 
                0.562866,0.165368,0.504692; 
                0.569172,0.167454,0.504105; 
                0.575490,0.169530,0.503466; 
                0.581819,0.171596,0.502777; 
                0.588158,0.173652,0.502035; 
                0.594508,0.175701,0.501241; 
                0.600868,0.177743,0.500394; 
                0.607238,0.179779,0.499492; 
                0.613617,0.181811,0.498536; 
                0.620005,0.183840,0.497524; 
                0.626401,0.185867,0.496456; 
                0.632805,0.187893,0.495332; 
                0.639216,0.189921,0.494150; 
                0.645633,0.191952,0.492910; 
                0.652056,0.193986,0.491611; 
                0.658483,0.196027,0.490253; 
                0.664915,0.198075,0.488836; 
                0.671349,0.200133,0.487358; 
                0.677786,0.202203,0.485819; 
                0.684224,0.204286,0.484219; 
                0.690661,0.206384,0.482558; 
                0.697098,0.208501,0.480835; 
                0.703532,0.210638,0.479049; 
                0.709962,0.212797,0.477201; 
                0.716387,0.214982,0.475290; 
                0.722805,0.217194,0.473316; 
                0.729216,0.219437,0.471279; 
                0.735616,0.221713,0.469180; 
                0.742004,0.224025,0.467018; 
                0.748378,0.226377,0.464794; 
                0.754737,0.228772,0.462509; 
                0.761077,0.231214,0.460162; 
                0.767398,0.233705,0.457755; 
                0.773695,0.236249,0.455289; 
                0.779968,0.238851,0.452765; 
                0.786212,0.241514,0.450184; 
                0.792427,0.244242,0.447543; 
                0.798608,0.247040,0.444848; 
                0.804752,0.249911,0.442102; 
                0.810855,0.252861,0.439305; 
                0.816914,0.255895,0.436461; 
                0.822926,0.259016,0.433573; 
                0.828886,0.262229,0.430644; 
                0.834791,0.265540,0.427671; 
                0.840636,0.268953,0.424666; 
                0.846416,0.272473,0.421631; 
                0.852126,0.276106,0.418573; 
                0.857763,0.279857,0.415496; 
                0.863320,0.283729,0.412403; 
                0.868793,0.287728,0.409303; 
                0.874176,0.291859,0.406205; 
                0.879464,0.296125,0.403118; 
                0.884651,0.300530,0.400047; 
                0.889731,0.305079,0.397002; 
                0.894700,0.309773,0.393995; 
                0.899552,0.314616,0.391037; 
                0.904281,0.319610,0.388137; 
                0.908884,0.324755,0.385308; 
                0.913354,0.330052,0.382563; 
                0.917689,0.335500,0.379915; 
                0.921884,0.341098,0.377376; 
                0.925937,0.346844,0.374959; 
                0.929845,0.352734,0.372677; 
                0.933606,0.358764,0.370541; 
                0.937221,0.364929,0.368567; 
                0.940687,0.371224,0.366762; 
                0.944006,0.377643,0.365136; 
                0.947180,0.384178,0.363701; 
                0.950210,0.390820,0.362468; 
                0.953099,0.397563,0.361438; 
                0.955849,0.404400,0.360619; 
                0.958464,0.411324,0.360014; 
                0.960949,0.418323,0.359630; 
                0.963310,0.425390,0.359469; 
                0.965549,0.432519,0.359529; 
                0.967671,0.439703,0.359810; 
                0.969680,0.446936,0.360311; 
                0.971582,0.454210,0.361030; 
                0.973381,0.461520,0.361965; 
                0.975082,0.468861,0.363111; 
                0.976690,0.476226,0.364466; 
                0.978210,0.483612,0.366025; 
                0.979645,0.491014,0.367783; 
                0.981000,0.498428,0.369734; 
                0.982279,0.505851,0.371874; 
                0.983485,0.513280,0.374198; 
                0.984622,0.520713,0.376698; 
                0.985693,0.528148,0.379371; 
                0.986700,0.535582,0.382210; 
                0.987646,0.543015,0.385210; 
                0.988533,0.550446,0.388365; 
                0.989363,0.557873,0.391671; 
                0.990138,0.565296,0.395122; 
                0.990871,0.572706,0.398714; 
                0.991558,0.580107,0.402441; 
                0.992196,0.587502,0.406299; 
                0.992785,0.594891,0.410283; 
                0.993326,0.602275,0.414390; 
                0.993834,0.609644,0.418613; 
                0.994309,0.616999,0.422950; 
                0.994738,0.624350,0.427397; 
                0.995122,0.631696,0.431951; 
                0.995480,0.639027,0.436607; 
                0.995810,0.646344,0.441361; 
                0.996096,0.653659,0.446213; 
                0.996341,0.660969,0.451160; 
                0.996580,0.668256,0.456192; 
                0.996775,0.675541,0.461314; 
                0.996925,0.682828,0.466526; 
                0.997077,0.690088,0.471811; 
                0.997186,0.697349,0.477182; 
                0.997254,0.704611,0.482635; 
                0.997325,0.711848,0.488154; 
                0.997351,0.719089,0.493755; 
                0.997351,0.726324,0.499428; 
                0.997341,0.733545,0.505167; 
                0.997285,0.740772,0.510983; 
                0.997228,0.747981,0.516859; 
                0.997138,0.755190,0.522806; 
                0.997019,0.762398,0.528821; 
                0.996898,0.769591,0.534892; 
                0.996727,0.776795,0.541039; 
                0.996571,0.783977,0.547233; 
                0.996369,0.791167,0.553499; 
                0.996162,0.798348,0.559820; 
                0.995932,0.805527,0.566202; 
                0.995680,0.812706,0.572645; 
                0.995424,0.819875,0.579140; 
                0.995131,0.827052,0.585701; 
                0.994851,0.834213,0.592307; 
                0.994524,0.841387,0.598983; 
                0.994222,0.848540,0.605696; 
                0.993866,0.855711,0.612482; 
                0.993545,0.862859,0.619299; 
                0.993170,0.870024,0.626189; 
                0.992831,0.877168,0.633109; 
                0.992440,0.884330,0.640099; 
                0.992089,0.891470,0.647116; 
                0.991688,0.898627,0.654202; 
                0.991332,0.905763,0.661309; 
                0.990930,0.912915,0.668481; 
                0.990570,0.920049,0.675675; 
                0.990175,0.927196,0.682926; 
                0.989815,0.934329,0.690198; 
                0.989434,0.941470,0.697519; 
                0.989077,0.948604,0.704863; 
                0.988717,0.955742,0.712242; 
                0.988367,0.962878,0.719649; 
                0.988033,0.970012,0.727077; 
                0.987691,0.977154,0.734536; 
                0.987387,0.984288,0.742002; 
                0.987053,0.991438,0.749504];
            num = size(raw,1);
            vec = linspace(0,num+1,N+2);
            map = interp1(1:num,raw,vec(2:end-1),'linear','extrap'); % Lab
            map = max(0,min(1,map));  
            
        end
        
    end
end
%------------------
function hm = getFaceHarmMean(G)
  tp =@(clambda) getTwoPointOperator(G);
  hmf =@(hflambda) getHarmonicAvgOpeartor(G);
  hm =@(clambda) hmf(tp(clambda));
end


%-------------------------------------------------------------------------%
function tp = getTwoPointOperator(G)
% Mappings from cells to its faces
cells = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
faces = G.cells.faces(:,1);
% Vector from cell to face centroid
C = G.faces.centroids(faces,:) - G.cells.centroids(cells,:);
% Oriented normals
sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1;
N   = bsxfun(@times, sgn, G.faces.normals(faces, :));
% Make function
cn  = sum(C.*N,2)./sum(C.*C,2);
tp = @(lambda) cn.*lambda(cells);
end

%-------------------------------------------------------------------------%
function ha = getHarmonicAvgOpeartor(G)
% Harmonig averaging operator
faces = G.cells.faces(:,1);
M = sparse(faces, 1:numel(faces), 1, G.faces.num, numel(faces));
ha = @(T) 1./(M*(1./T));
end
function allDiv = getAllDiv(G)
%% 
   nc = G.cells.num;
   nf = G.faces.num;
   Nall = G.faces.neighbors;
   internal = all(Nall>0,2);
   ifn = find(internal);
   efn = find(~internal);
   inf = numel(ifn);
   N=Nall(internal,:);
   Nb = Nall(~internal,:);
   signb = 2*(Nb(:,1)>0)-1;
   Nb=sum(Nb,2);
   C  = sparse([ifn; ifn], N, ones(inf, 1) * [1, -1], nf, nc);
   C  = C + sparse(efn, Nb, signb, nf, nc);
   allDiv =@(x) (C'*x)./G.cells.volumes;
end