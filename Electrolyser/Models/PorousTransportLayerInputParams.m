classdef PorousTransportLayerInputParams < ElectronicComponentInputParams
    
    properties
        
        solidVolumeFraction  % Solid volume fraction [-]
        leverettCoefficients  % Leverett coefficient that enters in the computation of the capillary pressure, see leverett.m 
        theta                % Water contact angle, enters in the computation of the capillary pressure [-]
        permeability         % Permeability [Darcy]
        tortuosity % Tortuosity
        species % species struct 
        % species.OH.molecularWeight    : Molecular weight of OH [kg mol^-1]
        % species.OH.V0    : Partial molar volume of OH [m^3 mol^-1]
        % species.OH.D     : Diffusion coefficient [m^2 s^-1]
        % species.OH.t     : Transference coefficient [-]
        % species.OH.chargeNumber     : Charge number [-]
        % species.K.molecularWeight     : Molecular weight of K [kg mol^-1]
        % species.K.V0     : Partial molar volume of K [m^3 mol^-1]
        % species.H2O.molecularWeight   : Molecular weight of H2O [kg mol^-1]
        % species.H2O.kLV  : Liquid-vapor exchange rate
        % species.H2O.mu0  : Standard chemical potential
        % species.H2O.V0   : Partial molar volume of H2O [m^3 mol^-1]

        externalCouplingTerm

        Boundary
        
    end
    
    methods
        
        function inputparams = PorousTransportLayerInputParams(jsonstruct)
            
            inputparams = inputparams@ElectronicComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
