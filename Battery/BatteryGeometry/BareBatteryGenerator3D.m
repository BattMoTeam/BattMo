classdef BareBatteryGenerator3D < BatteryGenerator
% setup 1D grid 
    properties
        
        sepnx = 10;
        nenx  = 10;
        penx  = 10;
        nenr  = 40; % discretization for solid diffusion
        penr  = 40; % discretization for solid diffusion
        fac   = 1;

        use_thermal
        
    end
    
    methods
        
        function gen = BareBatteryGenerator3D()
          gen = gen@BatteryGenerator();  
        end
            
        function paramobj = updateBatteryInputParams(gen, paramobj)
            if paramobj.include_current_collectors
                warning('This geometry does not include current collectors but input data has been given for those');
            end
            paramobj = gen.setupBatteryInputParams(paramobj, []);
        end
        
        function [paramobj, gen] = setupGrid(gen, paramobj, ~)
        % paramobj is instance of BatteryInputParams
        % setup paramobj.G
            sepnx = gen.sepnx;
            nenx  = gen.nenx;
            penx  = gen.penx;

            nxs = [nenx; sepnx; penx];

            % following dimension are those from Chen's paper
            xlength = 1e-5*[8.52; 1.2; 7.56];
            ylength = 1.58;
            zlength = 0.065;

            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            G = tensorGrid(x, [0; ylength], [0; zlength]);
            G = computeGeometry(G); 

            paramobj.G = G;
            gen.G = G;
            
        end

        function gen = applyResolutionFactors(gen)
            
            fac = gen.fac;
            
            gen.sepnx  = gen.sepnx*fac;
            gen.nenx   = gen.nenx*fac;
            gen.penx   = gen.penx*fac;
            
        end
            
        function paramobj = setupElectrolyte(gen, paramobj, params)
            
            params.cellind =  (1 : (gen.nenx + gen.sepnx + gen.penx))';
            params.Separator.cellind = gen.nenx + (1 : gen.sepnx)';
            
            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);
            
        end
        
        function paramobj = setupElectrodes(gen, paramobj, params)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            sd  = 'SolidDiffusion';
            
            sepnx = gen.sepnx; 
            nenx  = gen.nenx; 
            penx  = gen.penx; 
            
            %% parameters for negative electrode

            params.(ne).cellind = (1 : nenx)';
            params.(ne).(am).cellind = params.(ne).cellind;
            
            % boundary setup for negative electrode
            params.(ne).(am).bcfaces = 1;
            params.(ne).(am).bccells = 1;
            
            %% parameters for positive electrode
            
            pe_indstart = nenx + sepnx;
            params.(pe).cellind =  pe_indstart + (1 :  penx)';
            params.(pe).(am).cellind = params.(pe).cellind;
            
            % boundary setup for positive electode
            params.(pe).(am).bcfaces = penx + 1;
            params.(pe).(am).bccells = penx;

            paramobj = setupElectrodes@BatteryGenerator(gen, paramobj, params);

            % N and np are parameters for the full diffusion model
            if strcmp(paramobj.(ne).(am).diffusionModelType, 'full')
                paramobj.(ne).(am).(sd).N  = gen.nenr;
                paramobj.(ne).(am).(sd).np = gen.nenx;
            end
            
            if strcmp(paramobj.(pe).(am).diffusionModelType, 'full')
                paramobj.(pe).(am).(sd).N  = gen.penr;
                paramobj.(pe).(am).(sd).np = gen.penx;
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
