classdef SiliconGraphiteActiveMaterialInputParams < ElectronicComponentInputParams
%
% Input parameter class for :class:`ActiveMaterial <Electrochemistry.ActiveMaterial>`
% 
    properties

        Graphite

        Silicon

        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)
        
        volumeFraction % Volume fraction of the whole material (binder and so on included)
        

    end

    methods

        function paramobj = SiliconGraphiteActiveMaterialInputParams(jsonstruct)

            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);

            paramobj.Graphite = ActiveMaterialInputParams(pick('Graphite'));
            paramobj.Silicon = ActiveMaterialInputParams(pick('Silicon'));
            
            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)

            gr = 'Graphite';
            si = 'Silicon';

            paramobj = mergeParameters(paramobj,  {gr, 'volumeFraction'}, {'volumeFraction'});
            paramobj = mergeParameters(paramobj,  {si, 'volumeFraction'}, {'volumeFraction'});

            paramobj = validateInputParams@ElectronicComponentInputParams(paramobj);
            
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
