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

        function state = updateConductivity(model, state)
            
            cLi = state.cs{1};
            T = state.T;
            
            func = model.updateConductivityFunc;
            
            state.conductivity = func(cLi, T);
        end
        
        function state = updateChemicalCurrent(model, state)
            
            cLi          = state.cs{1}; % concentration of Li+
            T            = state.T;     % temperature
            phi          = state.phi;   % potential
            conductivity = state.conductivity;   % potential
            
            cs  = state.cs;         
            
            ncomp = model.ncomp; % number of components
            sp = model.sp;

            % calculate the concentration derivative of the chemical potential for each species in the electrolyte
            R = model.constants.R;
            for ind = 1 : ncomp
                dmudcs{ind} = R .* T ./ cs{ind};
            end
                        
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
            
            state.dmudcs = dmudcs;
            state.jchems = jchems;
            
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

