classdef GasDiffusionCellGridGenerator


    properties

        % Global grid
        % It is stored here because it is shared by many setup functions
        parentGrid

    end

    methods

        function [inputparams, gen] = updateGasDiffusionCellInputParams(gen, inputparams, params)
        % this function is the main class function as it returns an updated inputparams object with grid structure
            error('virtual function');
        end

        function [inputparams, gen] = setupGasDiffusionCellInputParams(gen, inputparams, params)
        % main function : add grid and coupling to inputparams structure

            [inputparams, gen] = gen.setupGrid(inputparams, params);
            inputparams = gen.setupExternalCoupling(inputparams, params);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
            error('virtual function');
        end
        

        function inputparams = setupExternalCoupling(gen, inputparams, params)

            ncontrols = params.nControls;

            for icontrol = 1 : ncontrols
                
                coupTerm = couplingTerm(sprintf('controlElement-%d', icontrol), {'boundary'});
                
                coupTerm.couplingfaces = params.bcfaces{icontrol};
                coupTerm.couplingcells = params.bccells{icontrol};

                inputparams.externalCouplingTerms{icontrol} = coupTerm;

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
