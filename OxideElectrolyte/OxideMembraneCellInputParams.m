classdef OxideMembraneCellInputParams < ComponentInputParams
    
    properties

        T
        
        Anode
        Cathode
        Electrolyte
        Control
        
        couplingTerms

        dx
        farea
        
    end
    
    methods
        
        function inputparams = OxideMembraneCellInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
    
            pick = @(fd) pickField(jsonstruct, fd);
            inputparams.(an)    = OxideMembraneElectrodeInputParams(pick(an));
            inputparams.(ct)    = OxideMembraneElectrodeInputParams(pick(ct));
            inputparams.(elyte) = OxideMembraneElectrolyteInputParams(pick(elyte));
            inputparams.(ctrl)  = OxideMembraneControlInputParams(pick(ctrl));

            inputparams = mergeParameters(inputparams, {{'T'}       , ...
                                                  {elyte, 'T'}, ...
                                                  {an, 'T'}   , ...
                                                  {ct, 'T'}});
            
            inputparams = mergeParameters(inputparams, {{elyte, 'muEl0'}, ...
                                                  {ct, 'muEl0'}   , ...
                                                  {an, 'muEl0'}});

            inputparams = mergeParameters(inputparams, {{elyte, 'Keh'}, ...
                                                  {ct, 'Keh'}   , ...
                                                  {an, 'Keh'}});


            inputparams.couplingTerms = {};
            
        end
        
    end
    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
