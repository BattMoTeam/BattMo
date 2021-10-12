classdef BinaryElectrolyte < Electrolyte
%
% Binary Electrolyte model
    
    properties

        ionicConductivityFittingCoefficients
        diffusionConcentrationFittingCoefficients
        diffusionTemperatureFittingCoefficients
        conductivityFactor = 1e-4
        
    end
    
    methods
        
        function model = BinaryElectrolyte(paramobj)
            
            model = model@Electrolyte(paramobj);
            fdnames = {'ionicConductivityFittingCoefficients', 
                       'diffusionConcentrationFittingCoefficients',
                       'diffusionTemperatureFittingCoefficients'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function state = updateConcentrations(model, state)
            state.cs{2} = state.cs{1};
        end
        
        function state = updateChemicalCurrent(model, state)
            
            cLi = state.cs{1}; % concentration of Li+
            T   = state.T;     % temperature
            phi = state.phi;   % potential
            
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
            cnst = model.ionicConductivityFittingCoefficients;            
           
            c = cLi;

            % Ionic conductivity, [S m^-1]
            conductivity = model.conductivityFactor.* c .*(polyval(cnst(end:-1:1, 1), c) + ...
                                                           polyval(cnst(end:-1:1, 2), c) .* T + ...
                                                           polyval(cnst(end:-1:1, 3), c) .* T.^2).^2;

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
            cnst = model.diffusionConcentrationFittingCoefficients;
            Tgi = model.diffusionTemperatureFittingCoefficients;
            
            % Diffusion coefficient, [m^2 s^-1]
            D = 1e-4 .* 10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ( T - Tgi(1) - Tgi(2) .* c .* 1e-3) + cnst(2, 1) .* ...
                                  c .* 1e-3) );
            
            % calculate the effective diffusion coefficients in porous media
            state.D = D .* model.volumeFraction .^1.5;
        
        end
        
    end
end

