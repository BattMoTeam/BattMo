classdef orgLiPF6 < Electrolyte
    
    properties

        conductivityFactor = 1e-4
        
    end
    
    methods
        
        function model = orgLiPF6(paramobj)
            
            model = model@Electrolyte(paramobj);
            fdnames = {'conductivityFactor'};

            model = dispatchParams(model, paramobj, fdnames);
        end
        
        function state = updateChemicalCurrent(model, state)
            
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
            state.dmudcs = dmudcs;
                        
            %% Calculate transport parameters
            % Calculate ionic conductivity consntants for the ionic conductivity calculation
            cnst = [-10.5   , 0.074    , -6.96e-5; ...
                    0.668e-3, -1.78e-5 , 2.80e-8; ...
                    0.494e-6, -8.86e-10, 0];            
            
            c = cLi;

            % Ionic conductivity, [S m^-1]
            conductivity = model.conductivityFactor.* c .*(...
                polyval(cnst(end:-1:1,1),c) + ...
                polyval(cnst(end:-1:1,2),c) .* T + ...
                polyval(cnst(end:-1:1,3),c) .* T.^2).^2;

            % volume fraction of electrolyte
            volfrac = model.volumeFraction;
            % Compute effective ionic conductivity in porous media
            conductivityeff = conductivity .* volfrac .^1.5;
            
            % setup chemical fluxes
            jchems = cell(1, ncomp);
            F = model.constants.F;
            for i = 1 : ncomp
                coeff = conductivityeff .* sp.t(i) .* dmudcs{i} ./ (sp.z(i).*F);
                jchems{i} = assembleFlux(model, cs{i}, coeff);
            end
            
            state.jchems = jchems;
            state.conductivity  = conductivityeff;
            
        end


        function state = updateDiffusionCoefficient(model, state)
            
            c = state.cs{1};
            T = state.T;
            
            % Calculate diffusion coefficients constant for the diffusion coefficient calcuation
            cnst = [ -4.43, -54;
                     -0.22, 0.0 ];

            Tgi = [ 229;
                    5.0 ];
            
            % Diffusion coefficient, [m^2 s^-1]
            D = 1e-4 .* 10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ( T - Tgi(1) - Tgi(2) .* c .* 1e-3) + cnst(2,1) .* ...
                                  c .* 1e-3) );
            
            % calculate the effective diffusion coefficients in porous media
            state.D = D .* model.volumeFraction .^1.5;
        
        end
        
    end
end

