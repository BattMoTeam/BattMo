classdef SolidDiffusionModel < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();

        %% Input parameters

        % Standard input parameters
        
        particleRadius                % the characteristic radius of the particle (symbol: rp)
        activationEnergyOfDiffusion   % the Arrhenius-type activation energy for diffusion (symbol: EaD)
        referenceDiffusionCoefficient % the pre-exponential reference diffusion coefficient in an Arrhenius-type equation (symbol: D0)
        volumetricSurfaceArea         % surface area of the active material - electrolyte interface per volume of electrode

        % Advanced input parameters
        volumeFraction % the ratio of the volume of the active material to the total volume
        
    end

    methods

        function model = SolidDiffusionModel(paramobj)

            model = model@BaseModel();

             % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);

            fdnames = {'particleRadius'               , ...
                       'activationEnergyOfDiffusion'  , ...
                       'referenceDiffusionCoefficient', ...
                       'volumetricSurfaceArea'        , ...
                       'volumeFraction'};

            model = dispatchParams(model, paramobj, fdnames);
        
        end
        
        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Temperature
            varnames{end + 1} = 'T';
            % Concentration at surface
            varnames{end + 1} = 'cSurface';
            % Diffusion coefficient
            varnames{end + 1} = 'D';
            % Volumetric reaction rate in mol/(s*m^3)
            varnames{end + 1} = 'Rvol';
            % Mass accumulation term
            varnames{end + 1} = 'massAccum';
            % Mass source term
            varnames{end + 1} = 'massSource';
            % Mass conservation equation
            varnames{end + 1} = 'massCons';

            model = model.registerVarNames(varnames);

            fn = @SolidDiffusionModel.updateDiffusionCoefficient;
            inputnames = {'T'};
            model = model.registerPropFunction({'D', fn, inputnames});

        end
        
        function state = updateDiffusionCoefficient(model, state)

            Tref = 298.15;  % [K]

            T = state.T;

            R   = model.constants.R;
            D0  = model.referenceDiffusionCoefficient;
            EaD = model.activationEnergyOfDiffusion;

            % Calculate solid diffusion coefficient, [m^2 s^-1]
            D = D0.*exp(-EaD./R*(1./T - 1/Tref));

            state.D = D;
            
        end
    
    end
    
end


%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
    

