classdef SolidDiffusionModelInputParams < InputParams 
%
% Base class for the solid diffusion models, see :class:`FullSolidDiffusionModelInputParams <Electrochemistry.FullSolidDiffusionModelInputParams>` and :class:`SimplifiedSolidDiffusionModelInputParams <Electrochemistry.SimplifiedSolidDiffusionModelInputParams>`
%
    properties

        %% Standard input parameters
        
        particleRadius                % the characteristic radius of the particle (symbol: rp)
        activationEnergyOfDiffusion   % the Arrhenius-type activation energy for diffusion (symbol: EaD)
        referenceDiffusionCoefficient % the pre-exponential reference diffusion coefficient in an Arrhenius-type equation (symbol: D0)
        volumetricSurfaceArea         % surface area of the active material - electrolyte interface per volume of electrode

        %% Advanced parameters
        
        volumeFraction % the ratio of the volume of the active material to the total volume (including porous space)
        
    end
    
    methods
        
        function inputparams = SolidDiffusionModelInputParams(jsonstruct)
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
