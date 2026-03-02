classdef LithiumPlatingLatzInputParams < InputParams
%
% Input parameter class for lithium platting. Reference:
% @article{hein2020electrochemical,
%   title     = {An electrochemical model of lithium plating and stripping in lithium ion batteries},
%   author    = {Hein, Simon and Danner, Timo and Latz, Arnulf},
%   journal   = {ACS Applied Energy Materials},
%   volume    = {3},
%   number    = {9},
%   pages     = {8519--8531},
%   year      = {2020},
%   publisher = {ACS Publications}
% }

    properties

        symmetryFactorPlating               % Symmetry factor for lithium plating/stripping reaction (typically ≈ 0.3)
        symmetryFactorStripping             % Symmetry factor for stripping (may be 1-symmetryFactorPlating)
        symmetryFactorChemicalIntercalation % Symmetry factor for charge-neutral chemical intercalation of plated lithium

        reactionRatePlating                 % Reaction rate constant for lithium plating (N⁰₀ Plating)
        reactionRateChemicalIntercalation   % Reaction rate constant for chemical intercalation of plated lithium
        reactionRateDirectIntercalation     % Reaction rate constant for direct lithium-ion intercalation into graphite

        thresholdParameter                  % Phenomenological parameter: minimum lithium amount needed to activate metal activity (see eqn (8))
        limitAmount                         % Limit amount of plated lithium corresponding to one monolayer on graphite surface (see eqn (26))
        platedReferenceConcentration        % Reference concentration for plated lithium [mol/m^3]

        platedLiMonolayerThickness          % Thickness of one monolayer of plated lithium around a particle [m]. See equation (S-3)

        volumetricSurfaceArea               % Interfacial surface area between graphite and electrolyte per unit volume [m²/m³]
        particleRadius                      % Radius of graphite particles, used in solid-state diffusion modeling
        volumeFraction                      % Porositywrite

        % Following parameters are not supported
        %
        % useSEI
        % SEIFraction
        % SEImolarMass       
        % SEIdensity     
        % SEIinitialThickness  
        % SEIconductivity   

    end

    methods

        function inputparams = LithiumPlatingLatzInputParams(jsonstruct)

            inputparams = inputparams@InputParams(jsonstruct);

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
