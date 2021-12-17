classdef LNMO < ActiveMaterial
    
    methods

        function model = LNMO(paramobj)
            
            model = model@ActiveMaterial(paramobj);
            model.n = 1;
            
        end
        
        function state = updateOCP(model, state)
        % Calculate the equilibrium open cirucuit potential of nmc111 according to the model used by Torchio et al [1].

        % Calculate the lithiation of the active material. This is a simplification for the initial code! The "real" value of
        % theta should be calculated using the surface concentration of Li and the maximum lithium concentration:
        %
            load('LNMO_OCV_data.mat');
            refT = 298.15;  % [K]

            T = state.T;
            c = state.cElectrode;
            
            theta = c ./ model.Li.cmax;
            % model.theta = (model.theta0 - model.theta100) .* model.SOC + model.theta100;
            
            % Calculate the open-circuit potential at the reference temperature for the given lithiation
                   
            
            
            
                    
            % Calculate the entropy change at the given lithiation
            
%             coeff1_dUdT = [0.199521039        , ...
%                            - 0.928373822      , ...
%                            + 1.364550689000003, ...
%                            - 0.611544893999998];
%             
%             coeff2_dUdT = [1                  , ...
%                            - 5.661479886999997, ...
%                            + 11.47636191      , ... 
%                            - 9.82431213599998 , ...
%                            + 3.048755063];
%             
%             dUdT = -1e-3.*polyval(coeff1_dUdT(end:-1:1),theta)./ polyval(coeff2_dUdT(end:-1:1),theta);
%             
%             % Calculate the open-circuit potential of the active material
%             OCP = refOCP + (T - refT) .* dUdT;
            if isa(theta, 'GenericAD')
                state.OCP = interp1(theta_ref, OCV_ref, theta.val);
            else
                state.OCP = interp1(theta_ref, OCV_ref, theta);
            end
            
        end
        
    end
end

%% References
%   [1] Versaci Daniele <daniele.versaci@polito.it> 
