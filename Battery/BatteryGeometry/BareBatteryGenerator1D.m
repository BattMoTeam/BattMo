classdef BareBatteryGenerator1D < BareBatteryGenerator
% setup 1D grid 
    properties
        
        sepnx  = 10;
        nenx   = 10;
        penx   = 10;
        fac = 1;
        
    end
    
    methods
        
        function gen = BareBatteryGenerator1D()
          gen = gen@BareBatteryGenerator();  
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

            nxs = [nenx; sepnx; penx];

            xlength = 1e-6*[100; 50; 80];
            ylength = 1e-2;

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
            
            params.cellind =  (1 : (gen.nenx + gen.sepnx + gen.penx))';
            params.Separator.cellind = gen.nenx + (1 : gen.sepnx)';
            
            paramobj = setupElectrolyte@BareBatteryGenerator(gen, paramobj, params);
        end
        
        function paramobj = setupElectrodes(gen, paramobj, params)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            
            sepnx = gen.sepnx; 
            nenx = gen.nenx; 
            penx = gen.penx; 
            
            %% parameters for negative electrode

            params.(ne).cellind = (1 : nenx)';
            
            % boundary setup for negative electrode
            params.(ne).bcfaces = 1;
            params.(ne).bccells = 1;
            
            %% parameters for positive electrode
            
            pe_indstart = nenx + sepnx;
            params.(pe).cellind =  pe_indstart + (1 :  penx)';
            
            % boundary setup for positive electode
            params.(pe).bcfaces = penx + 1;
            params.(pe).bccells = penx;
            
            paramobj = setupElectrodes@BareBatteryGenerator(gen, paramobj, params);

        end            I
                
    end
    
end





%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
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
