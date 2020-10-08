classdef orgLiPF6 < FvModel
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
        name    % Name of the electrolyte
        sp      % Dissolved species data structure
        ion     % Dissolved ionic species data structure
        solv    % Solvent data structure
        
        % Physical constants
        con = physicalConstants();
        
        % State properties
        c       % Molarity,                 [mol m^-3]
        m       % Molality,                 [mol kg^-1]
        wtp     % Weight percentace,        [wt%]
        T       % Temperature,              [K]
        eps     % Volume fraction,          [-]
        phi     % Electrostatic potential,  [V]
        IoSt    % Ionic strength
        j       % Ionic current density
        jchem   % Chemical component of ionic current density
        
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
        
        % Mesh properties
        X
        Xb
        N
        dombin
        
        
        % Finite volume solution properties
        chargeCont
        
        % number of components
        ncomp % 2 components : Li and PF6
    end
    
    methods
        function obj = orgLiPF6(c, T)
            %orgLiPF6 Construct an instance of the orgLIPF6 class
            %   obj = orgLiPF6(c, T) c is the concentration of LiPF6 with
            %   units [mol m^-3] and T is electrolyte temperature in Kelvin
            %   [K].
            obj.name = 'LiPF6';
            obj.c = c;
            obj.T = T;
            obj.ncomp = 2; 
            
            obj.sp.Li.c = c;
            obj.sp.PF6.c = c;
            
            % Set constant values
            obj.sp.Li.t = 0.399;           % Li+ transference number, [-]
            obj.sp.Li.z = 1;
            
            obj.sp.PF6.t = 1 - obj.sp.Li.t;   % PF6- transference number, [-]
            obj.sp.PF6.z = -1;
            
            obj.update();
        end
        
        function update(obj)
            obj.ionicQuantities();
            obj.conductivity();
            obj.diffusion();
        end
        
        function ionicQuantities(obj)
            
            obj.ion.cvec{1} = obj.sp.Li.c;
            obj.ion.cvec{2} = obj.sp.PF6.c;
            obj.ion.tvec{1} = obj.sp.Li.t .* ones(size(obj.sp.Li));
            obj.ion.tvec{2} = obj.sp.PF6.t .* ones(size(obj.sp.PF6));
            obj.ion.zvec{1} = obj.sp.Li.z;
            obj.ion.zvec{2} = obj.sp.PF6.z;
            for i = 1 : 2
                obj.ion.dmudc{i} = obj.con.R .* obj.T ./ obj.ion.cvec{i};
            end
            IoSt = 0.5 .* obj.ion.cvec{1}.*obj.ion.zvec{1}.^2./1000;
            IoSt = IoSt + 0.5 .* obj.ion.cvec{2}.*obj.ion.zvec{2}.^2./1000;
            obj.IoSt = IoSt;
            
        end
        
        function conductivity(obj, varargin)
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
                    obj.kappa = 1e-4 .* obj.c .* (...
                        (cnst(1,1) + cnst(2,1) .* obj.c + cnst(3,1) .* obj.c.^2) + ...
                        (cnst(1,2) + cnst(2,2) .* obj.c + cnst(3,2) .* obj.c.^2) .* obj.T + ...
                        (cnst(1,3) + cnst(2,3) .* obj.c) .* obj.T.^2) .^2;
            
        end
        
        function diffusion(obj)
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
                    obj.sp.Li.D = 1e-4 .* ...
                        10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ...
                                ( obj.T - Tgi(1) - Tgi(2) .* obj.c .* 1e-3) + ...
                                  cnst(2,1) .* obj.c .* 1e-3) );
            
        end
        
        
    end
end

%% References
%
%   [1] Journal ofThe Electrochemical Society, 152 (5) A882-A891 (2005),
%   DOI: 10.1149/1.1872737


