classdef orgLiPF6 < ComponentModel
%orgLiPF6 An electrolyte class for electrochemical modelling
%   The orgLiPF6 class describes the properties and
%   parameterization for organic electrolytes featuring lithium
%   hexafluorophosphate (LiPF6) salt dissolved in an alykl-carbonate
%   solvent. Common solvent materials are: 
%
%   PC          Propylene Carbonate     (ChemSpider ID: 7636)
%   EC          Ethylene Carbonate      (ChemSpider ID: 7030)
%   EMC         Ethyl Methyl Carbonate  (ChemSpider ID: 455390)
%   DMC         Dimethyl Carbonate      (ChemSpider ID: 11526)
%
%   The class calculates electrolyte properties based on experimental
%   parameterization studies described in the scientific literature.
%   The validity of the parameterization is limited to the conditions
%   in which it was reported.
%
%   Author: Simon Clark (simon.clark@sintef.no)
%   Usage:  This code is free to use "as-is" for the purpose of 
%           research at SINTEF without warranty of any kind. The code
%           is provided with the hope that it will be helpful. The 
%           author assumes no liability.         
%
%   Acknowledgement: This code builds on the work of many other
%   scientists over decades of research. Their work is gratefully
%   acknowledged and cited throughout the code. 
%
%   Revision History:
%       03.06.2020: SC (simon.clark@sintef.no) - New Energy Solutions 
%                   Initial version (0.0-alpha)
    
    properties
        % Identification properties
        sp      % Dissolved species data structure
        ion     % Dissolved ionic species data structure
        solv    % Solvent data structure

        % number of components
        
        compnames % 2 components in this implementation : Li and PF6 (can be generalized)
        ncomp % number of components
        
        
        % Physical constants
        con = physicalConstants();

        % Physicochemical properties
        eps         % Porosity
        rho         % Mass Density,                         [kg m^-3]
        mu          % Viscosity             
        kappaeff    % Porous media conductivity,            [S m^-1]
        lambda      % Thermal conductivity,                 [W m^-1 K^-1]
        lambdaeff   % Porous media thermal conductivity,    [W m^-1 K^-1]
        cp          % Heat Capacity
        sigma       % Surface Tension
        pvap        % Vapor Pressure    
        Deff        % Porous media diffusion coefficient,   [m^2 s^-1]
        
        % Finite volume solution properties
        chargeCont
        
    end
    
    methods
        
        function model = orgLiPF6(name, G, cells)
            
            % initialize as a ComponentModel in MRST
            model = model@ComponentModel(name);
            
            % generate the cell subgrid for the model
            model.G = genSubGrid(G, cells);
            
            % setup operators
            model.operators = localSetupOperators(model.G);
            
            % initialize the components of the model
            model.compnames = {'Li', 'PF6'};
            model.ncomp = numel(model.compnames);
            
            % define primary variables of the model
            pnames = {'phi', 'cLi'}; 
            model.pnames = pnames;
            
            % define state variables of the model
            names = {'phi'       , ...  % Electric potential
                     'T'         , ...  % Temperature
                     'cs'        , ...  % Concentrations
                     'LiSource'  , ...  % Li source term
                     'LiFlux'    , ...  % Li flux term
                     'chargeCont', ...  % Charge continuity
                     'massCont', ...    % Mass continuity
                    };
            model.names = names;
            
            model = model.setupVarDims();
            model.vardims('cs') = 2;

            % setup updating function for concentrations
            names = {'cs', 'chargeCont', 'LiFlux'};
            inputnames = {'cLi', 'T', 'LiSource'};
            updatefn = @(model, state) model.updateQuantities(state);
            
            for ind = 1 : numel(names)
                name = names{ind};
                if strcmp(name, 'cs')
                    name = VarName({'.'}, 'cs', 2);
                end
                model = model.addPropFunction(name, updatefn, inputnames, {'.'});
            end
            
            % add local aliases for each component
            aliases = {};
            for icn = 1 : model.ncomp
                compname = model.compnames{icn};
                name = sprintf('c%s', compname);
                lname = 'cs';
                varname = VarName({'.'}, lname, model.ncomp, icn);
                aliases{end + 1} = {name, varname};
            end
            
            model.aliases = aliases;

            % Set constant values
            [~, ind] = ismember('Li', model.compnames);
            tLi = 0.399;
            model.sp.t{ind} = tLi; % Li+ transference number, [-]
            model.sp.z{ind} = 1;
            
            [~, ind] = ismember('PF6', model.compnames);
            model.sp.t{ind} = 1 - tLi; % Li+ transference number, [-]
            model.sp.z{ind} = -1;
            
        end

        function state = initiateState(model, state)
            
            error('should not be used now');
            
            state = initiateState@ComponentModel(model, state);
            
            % instantiate the cell variables
            nc = model.G.cells.num;
            compnames = model.compnames;
            ncomp = model.ncomp;
            
            fieldnames = {'cs'};
            varnames = model.assignCurrentNameSpace(fieldnames);
            for ifn = 1 : numel(varnames)
                fieldname = varnames{ifn}.getfieldname;
                if isempty(state.(fieldname))
                    state.(fieldname) = cell(1, ncomp);
                end
            end
            
            
        end
        
        function state = updateQuantities(model, state)
            
     
            ncomp = model.ncomp;    % number of components
           
            cLi = state.cs{1};      % concentration of Li+
            T = state.T;            % temperature
            
            state.cs{2} = cLi;      % set counterion concentration?
            
            cs  = state.cs;         
            phi = state.phi;
            
            % calculate the concentration derivative of the chemical
            % potential for each species in the electrolyte
            for ind = 1 : ncomp
                dmudcs{ind} = model.con.R .* T ./ cs{ind};
            end
            
            % calculate ionic strength for a 2-component system
            IoSt = 0.5 .* cs{1}.*model.sp.z{1}.^2./1000;
            IoSt = IoSt + 0.5 .* cs{2}.*model.sp.z{2}.^2./1000;
            
            %% Calculate transport parameters
            % Calcualte ionic conductivity
            % consntants for the ionic conductivity calculation
            cnst = [-10.5   , 0.074    , -6.96e-5; ...
                    0.668e-3, -1.78e-5 , 2.80e-8; ...
                    0.494e-6, -8.86e-10, 0];            
            
            c = cLi;

            % Ionic conductivity, [S m^-1]
            kappa = 1e-4.* c .*(...
                polyval(cnst(end:-1:1,1),c) + ...
                polyval(cnst(end:-1:1,2),c) .* T + ...
                polyval(cnst(end:-1:1,3),c) .* T.^2).^2;
            
            % Calculate diffusion coefficients
            % constant for the diffusion coefficient calcuation
            cnst = [ -4.43, -54;
                     -0.22, 0.0 ];

            Tgi = [ 229;
                    5.0 ];
            
            % Diffusion coefficient, [m^2 s^-1]
            D = 1e-4 .* 10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ( T - Tgi(1) - Tgi(2) .* c .* 1e-3) + cnst(2,1) .* ...
                                  c .* 1e-3) );
            
            ncomp = model.ncomp;
            nf = model.G.faces.num;
            sp = model.sp;
            op = model.operators;
            
            % volume fraction of electrolyte
            eps = model.eps;
            
            % compute effective ionic conductivity in porous media
            kappaeff = kappa .* eps .^1.5;
            
            % setup chemical fluxes
            jchems = cell(1, ncomp);
            for i = 1 : ncomp
                coeff = kappaeff .* sp.t{i} .* dmudcs{i} ./ (sp.z{i}.*model.con.F);
                jchems{i} = op.harmFace(coeff).* op.Grad(cs{i});
            end
            
            jchem = jchems{1};
            for i = 2 : ncomp
                jchem = jchem + jchems{i};
            end
            
            % Ionic current density due to the electrochemical potential gradient
            j = op.harmFace(kappaeff).*(-1).*op.Grad(phi) - jchem;
            
            % We assume that LiSource has been setup
            LiSource = state.LiSource;
            
            % setup the charge continuity equation
            ind_Li = 1;
            chargeCont = - op.Div(j)./model.G.cells.volumes./model.con.F + LiSource.*model.sp.z{ind_Li};

            op = model.operators;
            
            % calculate the effective diffusion coefficients in porous
            % media
            Deff = D .* model.eps .^1.5;
            
            trans = op.harmFace(Deff);
            fluxDiff = - trans.*op.Grad(cLi);
            
            ind_Li = 1;
            fluxE = model.sp.t{ind_Li} ./ (model.sp.z{ind_Li} .* model.con.F) .* j;
            
            flux = fluxDiff + fluxE;
            
            state.LiFlux = flux;
            state.chargeCont = chargeCont;
            
        end

   
    end

    %% References
    %
    %   [1] Journal ofThe Electrochemical Society, 152 (5) A882-A891 (2005),
    %   DOI: 10.1149/1.1872737



end

