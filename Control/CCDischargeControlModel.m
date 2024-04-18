classdef CCDischargeControlModel < CCcontrolModel


    properties

        DRate
        
        rampupTime
        useCVswitch
        lowerCutoffVoltage

        % This value is initiated depending on DRate and battery model
        Imax

    end
    
    methods

        function model = CCDischargeControlModel(inputparams)
            
            model = model@CCcontrolModel(inputparams);
            
            fdnames = {'DRate'      , ...
                       'rampupTime' , ...
                       'useCVswitch', ...
                       'lowerCutoffVoltage'};

            model = dispatchParams(model, inputparams, fdnames);            

        end
        

        function func = setupStopFunction(model)

            func = setupStopFunction@ControlModel(model);
            
            if ~model.useCVswitch
                
                func = @(model, state, state_prev) (state.Control.E < model.Control.lowerCutoffVoltage);

            end
            
        end


        function  func = setupControlFunction(model)

            tup  = model.rampupTime;
            Umin = model.lowerCutoffVoltage;

            if model.useCVswitch
                
                func = @(time, I, E, Imax) rampupSwitchControl(time, ...
                                                               tup , ...
                                                               I   , ...
                                                               E   , ...
                                                               Imax, ...
                                                               Umin);
            else
                
                func = @(time, Imax) rampupControl(time, ...
                                                   tup , ...
                                                   Imax);
            end
                
        end
        
        function control = setupScheduleControl(model)

            control = setupScheduleControl@CCcontrolModel(model);
            control.CCDischarge = true;
            
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
