classdef SimplifiedSolidDiffusionModel < SolidDiffusionModel

    methods

        function model = SimplifiedSolidDiffusionModel(paramobj)

            model = model@SolidDiffusionModel(paramobj);

        end
        
        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@SolidDiffusionModel(model);

            varnames = {};
            % concentration
            varnames{end + 1} = 'cSurface';
            % concentration (over the whole particle)
            varnames{end + 1} = 'cAverage';
            % Solid diffusion equation
            varnames{end + 1} = 'solidDiffusionEq';
            
            model = model.registerVarNames(varnames);

            fn = @ActiveMaterial.updateDiffusionCoefficient;
            inputnames = {'T'};
            model = model.registerPropFunction({'D', fn, inputnames});

            fn = @ActiveMaterial.assembleSolidDiffusionEquation;
            inputnames = {'cSurface', 'cAverage', 'Rvol', 'D'};

            model = model.registerPropFunction({'solidDiffusionEq', fn, inputnames});
            
        end

        
        function state = assembleSolidDiffusionEquation(model, state)
        % We update the surface concentration of the charge carrier in the active material.
        % The surface concentration value is computed following polynomial method, as described in ref1 (see below)

            csurf = state.cSurface;
            caver = state.cAverage;
            D     = state.D;
            Rvol  = state.Rvol;

            rp  = model.rp;
            vsa = model.volumetricSurfaceArea;
            state.solidDiffusionEq = csurf - caver + (rp.*Rvol)./(5*vsa*D);
            
        end
    
    end
    
end


%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
    

