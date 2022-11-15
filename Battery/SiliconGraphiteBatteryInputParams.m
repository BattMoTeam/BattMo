classdef SiliconGraphiteBatteryInputParams < BatteryInputParams

    properties

        scenario = 'discharge'; % one of 'charge', 'first-charge', 'discharge'
        
    end
    
    methods
        
        function paramobj = SiliconGraphiteBatteryInputParams(jsonstruct)
            
            paramobj = paramobj@BatteryInputParams(jsonstruct);

            paramobj = paramobj.validateInputParams();
        end

        function paramobj = validateInputParams(paramobj)

            switch paramobj.scenario
                
              case 'first-charge'

                ne  = 'NegativeElectrode';
                pe  = 'PositiveElectrode';
                am  = 'ActiveMaterial';
                gr  = 'Graphite';
                si  = 'Silicon';
                itf = 'Interface';

                mats = {gr, si};

                for imat = 1 : numel(mats)
                    mat = mats{imat};
                    paramobj.(ne).(am).(mat).(itf).theta0 = 0.01;
                end
                
              case {'charge', 'discharge'}
                
                % do nothing
                
              otherwise
                
                error('scenario not recognized')
                
            end
                
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
