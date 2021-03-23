classdef orgLiPF6 < ElectrochemicalComponent
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

        % Physicochemical properties
        volumeFraction         % Porosity
        rho         % Mass Density,                         [kg m^-3]
        mu          % Viscosity             
        kappaeff    % Porous media conductivity,            [S m^-1]
        lambda      % Thermal conductivity,                 [W m^-1 K^-1]
        lambdaeff   % Porous media thermal conductivity,    [W m^-1 K^-1]
        cp          % Heat Capacity
        sigma       % Surface Tension
        pvap        % Vapor Pressure    
        Deff        % Porous media diffusion coefficient,   [m^2 s^-1]

        
        % names for book-keeping
        ionName
        ionFluxName 
        ionSourceName
        ionMassConsName
        ionAccumName
        
    end
    
    methods
        
        function model = orgLiPF6(G, cells)
            
            % initialize as a ComponentModel in MRST
            model = model@ElectrochemicalComponent();
            
            % generate the cell subgrid for the model
            model.G = genSubGrid(G, cells);
            
            % setup operators
            model.operators = localSetupOperators(model.G);
            
            % initialize the components of the model
            model.compnames = {'Li', 'PF6'};
            model.ncomp = numel(model.compnames);
            
            % define primary variables of the model

            % Set constant values
            [~, ind] = ismember('Li', model.compnames);
            tLi = 0.399;
            model.sp.t{ind} = tLi; % Li+ transference number, [-]
            model.sp.z{ind} = 1;
            
            [~, ind] = ismember('PF6', model.compnames);
            model.sp.t{ind} = 1 - tLi; % Li+ transference number, [-]
            model.sp.z{ind} = -1;
            
            % book-keeping variables
            model.ionName         = 'Li';
            model.ionFluxName     = 'LiFlux';
            model.ionSourceName   = 'LiSource';
            model.ionMassConsName = 'massCons';
            model.ionAccumName    = 'LiAccum';
            
        end


        function state  = updateCurrent(model, state) 
           
            cLi = state.cs{1}; % concentration of Li+
            T   = state.T;     % temperature
            phi = state.phi;   % potential
            
            state.cs{2} = cLi; % set counterion concentration?
            
            cs  = state.cs;         
            
            ncomp = model.ncomp; % number of components
            sp = model.sp;

            % calculate the concentration derivative of the chemical potential for each species in the electrolyte
            R = model.constants.R;
            for ind = 1 : ncomp
                dmudcs{ind} = R .* T ./ cs{ind};
            end
                        
            %% Calculate transport parameters
            % Calculate ionic conductivity consntants for the ionic conductivity calculation
            cnst = [-10.5   , 0.074    , -6.96e-5; ...
                    0.668e-3, -1.78e-5 , 2.80e-8; ...
                    0.494e-6, -8.86e-10, 0];            
            
            c = cLi;

            % Ionic conductivity, [S m^-1]
            kappa = 1e-4.* c .*(...
                polyval(cnst(end:-1:1,1),c) + ...
                polyval(cnst(end:-1:1,2),c) .* T + ...
                polyval(cnst(end:-1:1,3),c) .* T.^2).^2;

            % volume fraction of electrolyte
            eps = model.volumeFraction;
            % Compute effective ionic conductivity in porous media
            kappaeff = kappa .* eps .^1.5;
            
            % setup chemical fluxes
            jchems = cell(1, ncomp);
            F = model.constants.F;
            for i = 1 : ncomp
                coeff = kappaeff .* sp.t{i} .* dmudcs{i} ./ (sp.z{i}.*F);
                jchems{i} = assembleFlux(model, cs{i}, coeff);
            end
            
            j = assembleFlux(model, phi, kappaeff);
            for ind = 1 : ncomp
                j = j + jchems{ind};
            end

            state.j = j;
            
        end

        function state = updateIonFlux(model, state)
            state = model.updateLithiumFlux(state);
        end
        
        function state = updateLithiumFlux(model, state)
            
            % We assume that LiSource and current have been updated
            LiSource = state.LiSource;
            c = state.cs{1};
            j = state.j;
            T = state.T;
            
            %% 1 . Compute Flux from diffustion
            % Calculate diffusion coefficients constant for the diffusion coefficient calcuation
            cnst = [ -4.43, -54;
                     -0.22, 0.0 ];

            Tgi = [ 229;
                    5.0 ];
            
            % Diffusion coefficient, [m^2 s^-1]
            D = 1e-4 .* 10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ( T - Tgi(1) - Tgi(2) .* c .* 1e-3) + cnst(2,1) .* ...
                                  c .* 1e-3) );
            % calculate the effective diffusion coefficients in porous media
            Deff = D .* model.volumeFraction .^1.5;
            
            % Flux from diffusion
            fluxDiff = assembleFlux(model, c, Deff);
            
            %% 2. Compute Flux from electrical forces
            ind_Li = 1;
            F = model.constants.F;
            fluxE = model.sp.t{ind_Li} ./ (model.sp.z{ind_Li} .* F) .* j;
            
            %% 3. Sum the two flux contributions
            flux = fluxDiff + fluxE;
            state.LiFlux = flux;
           
        end
   
        function state = updateMassConservation(model, state)
            
            ionName         = model.ionName;
            ionFluxName     = model.ionFluxName;
            ionSourceName   = model.ionSourceName;
            ionAccumName    = model.ionAccumName;
            ionMassConsName = model.ionMassConsName;
            
            flux   = state.(ionFluxName);
            source = state.(ionSourceName);
            accum  = state.(ionAccumName);
            bcflux = 0;
            
            masscons = assembleConservationEquation(model, flux, bcflux, source, accum);
            
            state.(ionMassConsName) = masscons;
            
        end
    end

    %% References
    %
    %   [1] Journal ofThe Electrochemical Society, 152 (5) A882-A891 (2005),
    %   DOI: 10.1149/1.1872737



end

