classdef orgLiPF6 < SimpleModel
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
        rho         % Mass Density,                         [kg m^-3]
        mu          % Viscosity             
        kappa       % Conductivity,                         [S m^-1]
        kappaeff    % Porous media conductivity,            [S m^-1]
        lambda      % Thermal conductivity,                 [W m^-1 K^-1]
        lambdaeff   % Porous media thermal conductivity,    [W m^-1 K^-1]
        cp          % Heat Capacity
        sigma       % Surface Tension
        pvap        % Vapor Pressure    
        D           % Diffusion coefficient,                [m^2 s^-1]
        Deff        % Porous media diffusion coefficient,   [m^2 s^-1]
        
        % Finite volume solution properties
        chargeCont
        
    end
    
    methods
        
        function model = orgLiPF6(G, cells)
            model = model@SimpleModel();
            model.G = genSubGrid(G, cells);
            model.speciesnames = {'Li', 'PF6'};
            nodel.ncomp = numel(model.speciesnames);
        end
        
        function state = initializeState(model, state)

            varnames = model.getVarNames();
            for i = 1 : numel(varnames)
                varname = varnames{i};
                if ~isfield(state, varname)
                    state.(varname) = [];
                end
            end

            c = model.getProp(state, 'c-Li');
            state = model.setProp(state, 'c-PF6');
            
            % Set constant values
            [~, ind] = ismember('Li', compnames);
            tLi = 0.399;
            model.sp.t{ind} = tLi;           % Li+ transference number, [-]
            model.sp.z{ind} = 1;
            
            [~, ind] = ismember('PF6', compnames);
            model.sp.t{ind} = 1 - tLi;           % Li+ transference number, [-]
            model.sp.z{ind} = -1;
            
            model.update();
            
        end
        
        function [globalnames, localnames] = getModelPrimaryVarNames(model)
            localnames = {'phi', 'c-Li'}; % name 'c-Li' should match setup in getAffiliatedComponentNames
            globalnames = model.setupGlobalNames(localnames); 
        end
        
        function [globalnames, localnames] = getModelVarNames(model)
            
            [concnames, ionconcnames, jchemnames] = model.getAffiliatedComponentNames();

            localnames =  {m, ...    % Molality,                 [mol kg^-1]
                           wtp, ...  % Weight percentace,        [wt%]
                           eps, ...  % Volume fraction,          [-]
                           IoSt, ... % Ionic strength
                           j ...     % Ionic current density
                          };
            localnames = horzcat(concnames, ionconcnames, jchemnames, localnames);
            
            globalnames = model.setupGlobalNames(localnames);             
            
        end
        
        function [globalnames, localnames] = getVarNames(model)
        % this function 
            [globalnames, localnames] = model.getVarNames@SimpleModel();
            localnames                = horzcat(localnames, {'T'});
            globalnames               = horzcat(globalnames, {'T'});
        end
        
        
        function update(model)
            model.ionicQuantities();
            model.conductivity();
            model.diffusion();
        end
        
    end


    function [concnames, ionconcnames, jchemnames] = getAffiliatedComponentNames(model)
        compnames = model.compnames;
        concnames    = cellfun(@(x) sprintf('c-%s', x), compnames, 'uniformoutput', false);
        ionconcnames = cellfun(@(x) sprintf('ionc-%s', x), compnames, 'uniformoutput', false);            
        jchemnames   = cellfun(@(x) sprintf('jchem-%s', x), compnames, 'uniformoutput', false);
    end
    
    function state = ionicQuantities(model, state)
        
        ncomp = model.comp
        [concnames, ionconcnames, jchemnames] = getAffiliatedComponentNames(model)
        for i = 1 : numel(concnames)
            ionname = ionconcnames{i};
            cname = concnames{i};
            c = model.getProp(state, cname);
            state = model.setProp(state, ionname, c);
            tvec
            model.ion.tvec{1} = model.sp.Li.t .* ones(size(model.sp.Li));
            model.ion.tvec{2} = model.sp.PF6.t .* ones(size(model.sp.PF6));
            model.ion.zvec{1} = model.sp.Li.z;
            model.ion.zvec{2} = model.sp.PF6.z;
        end
        
        for i = 1 : 2
            model.ion.dmudc{i} = model.con.R .* model.T ./ model.ion.cvec{i};
        end
        IoSt = 0.5 .* model.ion.cvec{1}.*model.ion.zvec{1}.^2./1000;
        IoSt = IoSt + 0.5 .* model.ion.cvec{2}.*model.ion.zvec{2}.^2./1000;
        model.IoSt = IoSt;
        
    end
    
    function conductivity(model, varargin)
    %conductivity Calculates the ionic conductivity of the
    %eletrolyte in units [S m^-1].
    %   Electrolyte conductivity according to the model proposed by
    %   Val�en et al [1]. The model was made by performing a
    %   least-squares fit of experimental data with LiPF6 
    %   concenrations from 7.7e-6 M to 3.9 M and temperatures from
    %   263 K to 333 K. The solvent is 10 vol% PC, 27 vol% EC, 63
    %   vol% DMC.
        
    % Empirical fitting parameters
        cnst = [-10.5, 0.074, -6.96e-5; ...
                0.668e-3, -1.78e-5, 2.80e-8; ...
                0.494e-6, -8.86e-10, 0];
        
        % Electrolyte conductivity
        model.kappa = 1e-4 .* model.c .* (...
            (cnst(1,1) + cnst(2,1) .* model.c + cnst(3,1) .* model.c.^2) + ...
            (cnst(1,2) + cnst(2,2) .* model.c + cnst(3,2) .* model.c.^2) .* model.T + ...
            (cnst(1,3) + cnst(2,3) .* model.c) .* model.T.^2) .^2;
        
    end
    
    function diffusion(model)
    %diffusion Calculates the diffusion coefficient of Li+ ions in
    %the electrolyte in units [m2 s^-1].
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
        
        % Diffusion coefficient, [m^2 s^-1]
        model.sp.Li.D = 1e-4 .* ...
            10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ...
                      ( model.T - Tgi(1) - Tgi(2) .* model.c .* 1e-3) + ...
                      cnst(2,1) .* model.c .* 1e-3) );
        
    end
    
    
end

%% References
%
%   [1] Journal ofThe Electrochemical Society, 152 (5) A882-A891 (2005),
%   DOI: 10.1149/1.1872737


