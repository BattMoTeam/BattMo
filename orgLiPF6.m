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
            
            model = model@ComponentModel(name);
            model.G = genSubGrid(G, cells);
            
            % setup operators
            model.operators = localSetupOperators(model.G);
            
            model.compnames = {'Li', 'PF6'};
            model.ncomp = numel(model.compnames);
            
            % primary variables
            names = {'phi', 'c_Li'}; % name 'c_Li' should match the setup in getAffiliatedComponentNames
            model.pnames = names;
            
            % state variables
            names = {'phi', ...    % Potential
                     'T', ...      % Temperature
                     'cs', ...     % Concentrations
                     'ioncs', ...  % Ions concentrations (?)
                     'jchems', ... % Chemical fluxes
                     'dmudcs', ... % Chemical potentials (derivatives)
                     'm', ...      % Molality,              [mol kg^-1]
                     'kappa', ...  % Conductivity,          [S m^-1]
                     'LiSource', ...
                     'LiFlux', ...
                     'chargeCont', ...
                     'D', ...      % Diffusion coefficient, [m^2 s^-1]
                     'wtp', ...    % Weight percentace,     [wt%]
                     'IoSt', ...   % Ionic strength
                     'j' ...       % Ionic current density
                    };
            model.names = names;

            propfunctions = {};
            
            % setup updating function for concentrations
            name = 'cs';
            updatefn = @(model, state) model.updateConcentrations(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
            % setup updating function for dmudcs
            name = 'dmudcs';
            updatefn = @(model, state) model.updateIonicQuantities(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
            % setup updating function for ioncs            
            name = 'ioncs';
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
            % setup updating function for IoSt
            name = 'IoSt';
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
             % setup updating function for kappa
            name = 'kappa';
            updatefn = @(model, state) model.updateConductivity(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
            % setup updating function for D
            name = 'D';
            updatefn = @(model, state) model.updateDiffusion(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
            % setup updating function for jchem
            name = 'jchems';
            updatefn = @(model, state) model.updateChemicalFluxes(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;

            % setup updating function for j
            name = 'j';
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;

            % setup update property function for charge continuity (chargeCont)
            name = 'chargeCont';
            updatefn = @(model, state) model.updateChargeCont(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;

            % setup update property function for Li flux (LiFlux)
            name = 'LiFlux';
            updatefn = @(model, state) model.updateLiFlux(state);
            propfunction = PropFunction(name, updatefn, '.');
            propfunctions{end + 1} = propfunction;
            
            model.propfunctions = propfunctions;
            
            % add local aliases for each component
            aliases = {};
            fieldnames = {'c', 'ionc', 'jchem', 'dmudc'};
            for ifn = 1 : numel(fieldnames)
                fieldname = fieldnames{ifn};
                for icn = 1 : model.ncomp
                    compname = model.compnames{icn};
                    name = sprintf('%s_%s', fieldname, compname);
                    lname = sprintf('%ss', fieldname);
                    varname = VarName('.', lname, icn);
                    aliases{end + 1} = {name, varname};
                end
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

            state = initiateState@ComponentModel(model, state);
            
            % instantiate the cell variables
            nc = model.G.cells.num;
            compnames = model.compnames;
            ncomp = model.ncomp;
            
            fieldnames = {'cs', 'ioncs', 'jchems', 'dmudcs'};
            varnames = model.assignCurrentNameSpace(fieldnames);
            for ifn = 1 : numel(varnames)
                fieldname = varnames{ifn}.getfieldname;
                if isempty(state.(fieldname))
                    state.(fieldname) = cell(1, ncomp);
                end
            end
            
            
        end
        
        function state = initializeState(model, state)
            % used only in debugging for the moment
                
            state = model.initiateState(state);
            
            % instantiate cell variables
            nc = model.G.cells.num;
            compnames = model.compnames;
            ncomp = model.ncomp;
            
            fieldnames = {'cs', 'ioncs', 'jchems', 'dmudcs'};
            varnames = model.assignCurrentNameSpace(fieldnames);
            for ifn = 1 : numel(varnames)
                fieldname = varnames{ifn}.getfieldname;
                if isempty(state.(fieldname))
                    state.(fieldname) = cell(1, ncomp);
                end
            end
            
            [c, state] = model.getUpdatedProp(state, 'c_Li');

            state = model.setProp(state, 'c_PF6', c);

        end

        function state = updateConcentrations(model, state)
            [c_Li, state] = model.getUpdatedProp(state, 'c_Li');
            state = model.setProp(state, 'c_PF6', c_Li);
        end
        
        function state = updateIonicQuantities(model, state)
            
            ncomp = model.ncomp;
           
            [T, state]  = model.getUpdatedProp(state, 'T');
            [cs, state] = model.getUpdatedProp(state, 'cs');
            
            for ind = 1 : ncomp
                dmudcs{ind} = model.con.R .* T ./ cs{ind};
            end
            
            % this part is specific to 2 component system
            IoSt = 0.5 .* cs{1}.*model.sp.z{1}.^2./1000;
            IoSt = IoSt + 0.5 .* cs{2}.*model.sp.z{2}.^2./1000;
            
            state = model.setProp(state, 'dmudcs', dmudcs);
            state = model.setProp(state, 'ioncs', cs);
            state = model.setProp(state, 'IoSt', IoSt);
            
        end
        
        function state = updateConductivity(model, state)
        %   conductivity Calculates the ionic conductivity of the
        %   eletrolyte in units [S m^-1].
        %   Electrolyte conductivity according to the model proposed by
        %   Val�en et al [1]. The model was made by performing a
        %   least-squares fit of experimental data with LiPF6 
        %   concenrations from 7.7e-6 M to 3.9 M and temperatures from
        %   263 K to 333 K. The solvent is 10 vol% PC, 27 vol% EC, 63
        %   vol% DMC.
            
        % Empirical fitting parameters
            cnst = [-10.5   ,    0.074    ,    -6.96e-5; ...
                    0.668e-3,    -1.78e-5 ,    2.80e-8; ...
                    0.494e-6,    -8.86e-10,    0];
            
            % Electrolyte conductivity
            [T, state] = model.getUpdatedProp(state, 'T');
            [c, state] = model.getUpdatedProp(state, 'c_Li');
            
            kappa = 1e-4 .* c .* (...
                (cnst(1,1) + cnst(2,1) .* c + cnst(3,1) .* c.^2) + ...
                (cnst(1,2) + cnst(2,2) .* c + cnst(3,2) .* c.^2) .* T + ...
                (cnst(1,3) + cnst(2,3) .* c) .* T.^2) .^2;
            
            state = model.setProp(state, 'kappa', kappa);
            
        end
        
        function state = updateDiffusion(model, state)
        %   diffusion Calculates the diffusion coefficient of Li+ ions in
        %   the electrolyte in units [m2 s^-1].
        %   Diffusion coefficient according to the model proposed by
        %   Val�en et al [1]. The model was made by performing a
        %   least-squares fit of experimental data with LiPF6 
        %   concenrations from 7.7e-6 M to 3.9 M and temperatures from
        %   263 K to 333 K. The solvent is 10 vol% PC, 27 vol% EC, 63
        %   vol% DMC.
            
        % Empirical fitting parameters [1]
            cnst = [ -4.43, -54;
                     -0.22, 0.0 ];
            Tgi = [ 229;
                    5.0 ];
            
            [T, state] = model.getUpdatedProp(state, 'T');
            [c, state] = model.getUpdatedProp(state, 'c_Li');
            
            % Diffusion coefficient, [m^2 s^-1]
            D = 1e-4 .* 10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ( T - Tgi(1) - Tgi(2) .* c .* 1e-3) + cnst(2,1) .* ...
                                  c .* 1e-3) );
            
            state = model.setProp(state, 'D', D);
            
        end

        function state = updateChargeCont(model, state)
            
            op = model.operators;
            
            [j, state] = model.getUpdatedProp(state, 'j');
            [LiSource, state] = model.getUpdatedProp(state, 'LiSource');
            
            ind_Li = 1;
            chargeCont = - op.Div(j)./model.G.cells.volumes./model.con.F + LiSource.*model.sp.z{ind_Li};
                    
            state = model.setProp(state, 'chargeCont', chargeCont);
            
        end
        
        
        function state = updateLiFlux(model, state)
            
            op = model.operators;

            [D, state]   = model.getUpdatedProp(state, 'D');
            [cLi, state] = model.getUpdatedProp(state, 'c_Li');
            [j, state]   = model.getUpdatedProp(state, 'j');
            
            Deff = D .* model.eps .^1.5;
            
            trans = op.harmFace(Deff);
            fluxDiff = - trans.*op.Grad(cLi);
            
            ind_Li = 1;
            fluxE = model.sp.t{ind_Li} ./ (model.sp.z{ind_Li} .* model.con.F) .* j;
            
            flux = fluxDiff + fluxE;
            
            state = model.setProp(state, 'LiFlux', flux);
            
        end
        
        
        function state = updateChemicalFluxes(model, state);
            ncomp = model.ncomp;
            nf = model.G.faces.num;
            sp = model.sp;
            op = model.operators;
            
            [dmudcs, state] = model.getUpdatedProp(state, 'dmudcs');
            [cs, state]     = model.getUpdatedProp(state, 'cs');
            [phi, state]    = model.getUpdatedProp(state, 'phi');            
            [kappa, state]  = model.getUpdatedProp(state, 'kappa');

            eps = model.eps;
            
            % compute kappaeff
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
            
            state = model.setProp(state, 'jchems', jchems);
            state = model.setProp(state, 'j', j);
            
        end
        
   
    end

    %% References
    %
    %   [1] Journal ofThe Electrochemical Society, 152 (5) A882-A891 (2005),
    %   DOI: 10.1149/1.1872737



end

