classdef SwellingMaterial < ActiveMaterial
    
    properties

    end
    
    methods
        
        function model = SwellingMaterial(paramobj)
        %
        % ``paramobj`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %
            model = model@ActiveMaterial(paramobj);
            
            
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ActiveMaterial(model);
            
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            
            varnames = {'volumeFraction',
                        'radius',
                       {sd, 'cAverage'}};
            model = model.registerVarNames(varnames);

            fn = @SwellingMaterial.updateVolumeFraction;
            model = model.registerPropFunction({'volumeFraction', fn, {{sd, 'cAverage'}, 'radius'}});

            fn = @SwellingMaterial.SolidDiffusionModel.updateAverageConcentration;
            %% TODO check definition of update function for cAverage in FullSolidDiffusionModel (updateAverageConcentration)
            % If it can be used, then Xavier fix that
            model = model.registerPropFunction({{sd, 'cAverage'}, fn, {{sd, 'c'}}});
            
            fn = @SwellingMaterial.updateMassSource;
            model = model.registerPropFunction({{sd, 'massSource'}, fn, {{sd, 'Rvol'}, 'radius'}});

            fn = @SwellingMaterial.updateRadius;
            model = model.registerPropFunction({{'radius'}, fn, {{sd, 'cAverage'}}});


            %% TODO
            % update the functions that used volumefraction as a parameter before
            
        end

        function state = updateVolumeFraction(model, state)

            sd  = 'SolidDiffusion';

            c = state.(sd).cAverage;
            r = state.radius

            state.volumeFraction = ;
            
        end
        
        function state = updateMassCons(model, state)

            sd  = 'SolidDiffusion';
            
            Rvol = state.(sd).Rvol;
            r    = state.radius

            state.(sd).massSource = ;
            
        end
        
        function state = updateRadius(model, state)
            
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
