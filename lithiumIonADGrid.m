classdef lithiumIonADGrid < handle
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
        use_org
        AutoDiffBackend
    end
    
    methods
        function obj = lithiumIonADGrid(varargin)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin == 0
                SOC = 0.5;
                T   = 298.15;
            else
                SOC = varargin(1);
                T   = varargin(2);
            end
            obj.use_org=false;
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
            %obj.AutoDiffBackend = DiagonalAutoDiffBackend();
            %obj.AutoDiffBackend.useMex=false;
        end
        
        function spm(obj)
            %SPM Creates a single particle model of the Li-ion battery
            %   Detailed explanation goes here
            
            
        end

        function [dfdy, dfdyp] = odederfun(obj, t, y, yp) 
            res = obj.odefun(t, y, yp, 'useAD', true);
            if(numel(res.jac)==2)
                dfdy = res.jac{1};
                dfdyp = res.jac{2};
            else
                dfdy = res.jac{1}(:,1:res.numVars(1));
                dfdyp = res.jac{1}(:,res.numVars(1)+1:end);
            end
        end
    
         
        function [t, y] = p2d(obj)
           
            % Generate the FV mesh
            names           = {'NE';'SEP';'PE'};
            sizes           = [ 1e-6; 1e-6; 1e-6 ]*10;
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
                obj.(dd).operators.allDiv = getAllDiv(obj.(dd).Grid);
                %a =getTwoPointOperator(obj.(dd).Grid)
                %b=getHarmonicAvgOpeartor(obj.(dd).Grid)
                %obj.(dd).operators.harmFace =@(cvalue) b(a(cvalue));
                obj.(dd).operators.harmFace = getFaceHarmMean(obj.(dd).Grid);
                obj.(dd).operators.harmFaceBC =@(cvalue,faces) getFaceHarmBC(obj.(dd).Grid,cvalue,faces);
            end
            % Set initial conditions
            obj.icp2d();
            %% Is this constant?obj.elyte.sp.Li.Deff

            
            obj.elyte.sp.Li.Trans =  obj.elyte.operators.harmFace(obj.elyte.sp.Li.Deff);
            obj.ne.am.Li.Trans = obj.ne.operators.harmFace(obj.ne.am.Li.Deff);
            obj.pe.am.Li.Trans = obj.pe.operators.harmFace(obj.pe.am.Li.Deff);
            
            % Perform manual mesh chec
            if strcmpi(obj.check, 'on')==1
                figureWindow(obj)
                plotMesh(obj)
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
                    plotSummary(obj,t,y);
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
                coeff =  obj.elyte.kappaeff .* obj.elyte.ion.tvec{i} .* obj.elyte.ion.dmudc{i} ./ ...
                                     (obj.elyte.ion.zvec{i}.*obj.con.F);                                             
                jchems{i} = obj.elyte.operators.harmFace(coeff).* obj.elyte.operators.Grad(obj.elyte.ion.cvec{i});
                
            end
            obj.elyte.jchem = jchems{1};
            for i = 2 : ncomp
                obj.elyte.jchem = obj.elyte.jchem + jchems{i};
            
            end
            %   Ionic current density due to the electrochemical potential gradient
            obj.elyte.j = obj.elyte.operators.harmFace(obj.elyte.kappaeff).*(-1).*obj.elyte.operators.Grad(obj.elyte.phi) - obj.elyte.jchem;
      
            %% Electric current density
            % Active material NE

            obj.ne.j =  obj.ne.operators.harmFace(obj.ne.sigmaeff) ...
                            .* (-1) .* obj.ne.operators.Grad(obj.ne.am.phi);
               
               
            % add bc
             faceleft=1;
               bcLeft = 0;
               [t,cells] = obj.ne.operators.harmFaceBC(obj.ne.sigmaeff,faceleft)
               obj.ne.j_bcsource = obj.ne.am.phi*0.0;%NB hack to initialize zero ad
               obj.ne.j_bcsource(cells) =  t.*(bcLeft-obj.ne.am.phi(cells));
                                    
                % Active material PE                 
               obj.pe.j = obj.pe.operators.harmFace(obj.pe.sigmaeff)...
                            .* (-1) .* obj.pe.operators.Grad(obj.pe.am.phi);
               %ad bc
               faceright=obj.pe.Grid.faces.num;
               bcRight = obj.pe.E;
               [t,cells] = obj.pe.operators.harmFaceBC(obj.pe.sigmaeff,faceright);
               obj.pe.j_bcsource = obj.pe.am.phi*0.0;%NB hack to initialize zero ad
               obj.pe.j_bcsource(cells) =  t.*(bcRight-obj.pe.am.phi(cells));                       
            
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
            from = 1:obj.ne.N;
            to = 1:obj.ne.N;% not it is all
            elyte_ne=struct('from',from,'to',to);
            from = (obj.ne.N + obj.sep.N + 1):obj.elyte.N;
            to = 1:obj.pe.N;% now it is all
            elyte_pe=struct('from',from,'to',to);
            
            
            obj.ne.reactBV(obj.elyte.phi(elyte_ne.from));
            obj.pe.reactBV(obj.elyte.phi(elyte_pe.from));
            
            % Set source terms            
            %%%%% Li+ Sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Electrolyte NE Li+ source
            obj.elyte.sp.Li.source(elyte_ne.from) = +1 .* obj.ne.R;

            % Electrolyte PE Li+ source
            obj.elyte.sp.Li.source(elyte_pe.from) = -1 .* obj.pe.R;
            
            %%%%% Li0 Sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Active Material NE Li0 source
            obj.ne.am.Li.source(elyte_ne.to) = -1 .* obj.ne.R;

            % Active Material PE Li0 source
            obj.pe.am.Li.source(elyte_pe.to) = +1 .* obj.pe.R;
            
            %%%%% e- Sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Active Material NE Li0 source
            obj.ne.am.e.source(elyte_ne.to) = +1 .* obj.ne.R;

            % Active Material PE Li0 source
            obj.pe.am.e.source(elyte_pe.to) = -1 .* obj.pe.R;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Diffusion Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of diffusion mass flux
            %   Electrolyte Li+ Diffusion
            %%
            x = obj.elyte.sp.Li.c;
            flux = - obj.elyte.sp.Li.Trans.*obj.elyte.operators.Grad(x);
            obj.elyte.sp.Li.divDiff=  obj.elyte.operators.Div(flux)./obj.elyte.Grid.cells.volumes;
            
            %%
            x= obj.ne.am.Li.cs;
            flux = - obj.ne.am.Li.Trans.*obj.ne.operators.Grad(x);
            obj.ne.am.Li.divDiff=  obj.ne.operators.Div(flux)./obj.ne.Grid.cells.volumes;
            
            %%
            x = obj.pe.am.Li.cs;
            flux = - obj.pe.am.Li.Trans.*obj.pe.operators.Grad(x);
            obj.pe.am.Li.divDiff=  obj.pe.operators.Div(flux)./obj.pe.Grid.cells.volumes;
                
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Migration Flux                                           %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Divergence of the migration mass flux
            %   Electrolyte Li+ Migration
            flux = obj.elyte.sp.Li.t ./ (obj.elyte.sp.Li.z .* obj.con.F) ...
                                                .* obj.elyte.j;
                                
            obj.elyte.sp.Li.divMig   =  obj.elyte.operators.Div(flux)./obj.elyte.Grid.cells.volumes;                           

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

            obj.elyte.chargeCont = -(obj.elyte.operators.Div(  obj.elyte.j)./obj.elyte.Grid.cells.volumes) ./ obj.con.F+ ...
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
                                        
            %% Active material charge conobj.ne.operators.harmFaceBc(obj.ne.sigmaeff,faceleft);          

           obj.ne.am.e.chargeCont = ((-(obj.ne.operators.Div(  obj.ne.j)-obj.ne.j_bcsource))./obj.ne.Grid.cells.volumes)./-obj.con.F + ...
                                     obj.ne.am.e.source .* -1;
                                 
           %obj.pe.am.e.chargeCont = -div(  obj.pe.j, obj.pe.Xb) ./ -obj.con.F + ...
           %                          obj.pe.am.e.source .* -1;   
           obj.pe.am.e.chargeCont = ((-(obj.pe.operators.Div(obj.pe.j) - obj.pe.j_bcsource))./obj.pe.Grid.cells.volumes) ./ -obj.con.F + ...
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
       
        
    end
   
end
%------------------
function [T,cells] = getFaceHarmBC(G,cvalue,faces)
        cells=sum(G.faces.neighbors(faces,:));
        cn = sqrt(sum((G.faces.centroids(faces,:)- G.cells.centroids(cells,:)).^2,2));
        t = G.faces.areas(faces)./cn;
        T = t.*cvalue(cells);
end

function hm = getFaceHarmMean(G)
    internal = all(G.faces.neighbors>0,2);
    N = G.faces.neighbors(internal,:);
    ni = sum(internal);
    cd = sqrt(sum((G.cells.centroids(N(:,1),:) - G.cells.centroids(N(:,2),:)).^2,2));%NB
    t= G.faces.areas(internal)./cd
    A = sparse([[1:ni]';[1:ni]'],N,1,ni,G.cells.num);
    hm = @(cellvalue) 2.*t./(A*(1./cellvalue));
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
    
