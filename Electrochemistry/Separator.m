classdef Separator < BaseModel

    properties

        %% Input parameters

        % Standard parameters
        porosity             % the ratio of the volume free space to the total volume (symbol: varepsilon)
        density              % the mass density of the material (symbol: rho)
        bruggemanCoefficient % coefficient to determine effective transport parameters in porous media (symbol: beta)

        thermalConductivity  % Intrinsic Thermal conductivity of the electrolyte
        specificHeatCapacity % Specific Heat capacity of the electrolyte

        % Advanced parameters

        effectiveThermalConductivity
        effectiveVolumetricHeatCapacity

        %% Helper properties

        use_thermal

    end

    methods

        function model = Separator(inputparams)

            model = model@BaseModel();

            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            % in the case of the separator, probably this does not matter as no computation is actually done on this grid
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', true);

            fdnames = {'G'                              , ...
                       'porosity'                       , ...
                       'density'                        , ...
                       'bruggemanCoefficient'           , ...
                       'thermalConductivity'            , ...
                       'specificHeatCapacity'           , ...
                       'effectiveThermalConductivity'   , ...
                       'effectiveVolumetricHeatCapacity', ...
                       'use_thermal'};
            model = dispatchParams(model, inputparams, fdnames);

            if model.use_thermal

                if isempty(model.effectiveThermalConductivity)

                    bg = model.bruggemanCoefficient;
                    vf = 1 - model.porosity;

                    model.effectiveThermalConductivity = vf.^bg.*model.thermalConductivity;

                end

                if isempty(model.effectiveVolumetricHeatCapacity)

                    vf = 1 - model.porosity;

                    model.effectiveVolumetricHeatCapacity = vf.*model.density.*model.specificHeatCapacity;

                end

            end

        end

        function model = setTPFVgeometry(model, tPFVgeometry)
        % tPFVgeometry should be instance of TwoPointFiniteVolumeGeometry

            model.G.parentGrid.tPFVgeometry = tPFVgeometry;

        end

        function jsonstruct = exportParams(model)

            jsonstruct = exportParams@BaseModel(model);
            
            fdnames = {'porosity'                       , ...              
                       'density'                        , ...               
                       'bruggemanCoefficient'           , ...  
                       'thermalConductivity'            , ...   
                       'specificHeatCapacity'           , ...  
                       'effectiveThermalConductivity'   , ... 
                       'effectiveVolumetricHeatCapacity', ... 
                       'use_thermal'};
                       
            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                jsonstruct.(fdname) = model.(fdname);
            end

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
