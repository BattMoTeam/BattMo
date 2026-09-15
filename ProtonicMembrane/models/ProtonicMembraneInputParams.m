classdef ProtonicMembraneInputParams < ComponentInputParams
    
    properties

        T
        
        Anode
        Cathode
        Electrolyte
        Control
        
        couplingTerms

    end
    
    methods
        
        function inputparams = ProtonicMembraneInputParams(jsonstruct)

            jsonstruct = setDefaultStructField(jsonstruct, {'TimeStepping', 'useSwitch'}, true);
            jsonstruct = setDefaultStructField(jsonstruct, {'TimeStepping', 'fractionSwitch'}, 0.5);
            jsonstruct = setDefaultStructField(jsonstruct, {'TimeStepping', 'orderSwitch'}, 'I-first');
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
    
            pick = @(fd) pickField(jsonstruct, fd);
            inputparams.(an)    = ProtonicMembraneAnodeInputParams(pick(an));
            inputparams.(ct)    = ProtonicMembraneCathodeInputParams(pick(ct));
            inputparams.(elyte) = ProtonicMembraneElectrolyteInputParams(pick(elyte));
            inputparams.(ctrl)  = ProtonicMembraneControlInputParams(pick(ctrl));

            inputparams = mergeParameters(inputparams, {{'T'}       , ...
                                                  {elyte, 'T'}, ...
                                                  {an, 'T'}   , ...
                                                  {ct, 'T'}});
            
            inputparams = mergeParameters(inputparams, {{elyte, 'Ptot'}, ...
                                                  {an, 'Ptot'}   , ...
                                                  {ct, 'Ptot'}});

            inputparams = mergeParameters(inputparams, {{elyte, 'SU'}, ...
                                                  {an, 'SU'}});

            inputparams = mergeParameters(inputparams, {{elyte, 'steam_ratio'}, ...
                                                  {an, 'steam_ratio'}});

            inputparams = mergeParameters(inputparams, {{elyte, 'E_0'}, ...
                                                  {an, 'E_0'}});
            
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
