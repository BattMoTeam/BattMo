classdef BatteryGenerator1D < BatteryGenerator
% setup 1D grid 
    properties
        
        sepnx  = 10;
        nenx   = 10;
        penx   = 10;
        nenr   = 10; % discretization for solid diffusion
        penr   = 10; % discretization for solid diffusion
        ccnenx = 10;
        ccpenx = 10;
        fac = 1;
        include_current_collectors
        
    end
    
    methods
        
        function gen = BatteryGenerator1D()
          gen = gen@BatteryGenerator();  
        end
            
        function paramobj = updateBatteryInputParams(gen, paramobj)
            paramobj = gen.setupBatteryInputParams(paramobj, []);
        end
        
        function [paramobj, gen] = setupGrid(gen, paramobj, ~)
        % paramobj is instance of BatteryInputParams
        % setup paramobj.G
            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;
            ccnenx = gen.ccnenx;
            ccpenx = gen.ccpenx;

            nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
        
            %xlength = 1e-6*[10; 100; 50; 80; 10];
            afac=1;
            xlength = 1e-6*[25, 64*afac, 15, 57*afac, 15]';
    
            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            G = tensorGrid(x);
            G = computeGeometry(G); 

            paramobj.G = G;
            gen.G = G;
            
        end

        function gen = applyResolutionFactors(gen)
            
            fac = gen.fac;
            
            gen.sepnx  = gen.sepnx*fac;
            gen.nenx   = gen.nenx*fac;
            gen.penx   = gen.penx*fac;
            gen.ccnenx = gen.ccnenx*fac;
            gen.ccpenx = gen.ccpenx*fac;
            
        end
            
        function paramobj = setupElectrolyte(gen, paramobj, params)
            
            params.cellind = gen.ccnenx + (1 : (gen.nenx + gen.sepnx + gen.penx))';
            params.Separator.cellind = gen.ccnenx + gen.nenx + (1 : gen.sepnx)';
            
            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);
        end
        
        function paramobj = setupElectrodes(gen, paramobj, params)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            am = 'ActiveMaterial';
            
            sepnx = gen.sepnx; 
            nenx = gen.nenx; 
            penx = gen.penx; 
            ccnenx = gen.ccnenx; 
            ccpenx = gen.ccpenx;     
            
            %% parameters for negative electrode

            params.(ne).cellind = (1 : ccnenx + nenx)';
            params.(ne).(am).cellind = ccnenx + (1 : nenx)';
            params.(ne).(cc).cellind = (1 : ccnenx)';
            
            % boundary setup for negative current collector
            params.(ne).(cc).bcfaces = 1;
            params.(ne).(cc).bccells = 1;
            
            %% parameters for positive electrode
            
            pe_indstart = ccnenx + nenx + sepnx;
            params.(pe).cellind =  pe_indstart + (1 : ccpenx + penx)';
            params.(pe).(am).cellind = pe_indstart + (1 : penx)';
            params.(pe).(cc).cellind = pe_indstart + penx + (1 : ccpenx)';
            
            % boundary setup for positive current collector
            params.(pe).(cc).bcfaces = ccpenx + 1;
            params.(pe).(cc).bccells = ccpenx;
            
            paramobj = setupElectrodes@BatteryGenerator(gen, paramobj, params);

        end
                
        function paramobj = setupThermalModel(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
            params.couplingfaces = [];
            params.couplingcells = (1 : gen.G.cells.num)';
            paramobj = setupThermalModel@BatteryGenerator(gen, paramobj, params);
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
