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
        use_thermal
        
    end
    
    methods
        
        function gen = BatteryGenerator1D()
            gen = gen@BatteryGenerator();  
        end
            
        function paramobj = updateBatteryInputParams(gen, paramobj)
            gen.include_current_collectors = paramobj.include_current_collectors;
            gen.use_thermal = paramobj.use_thermal;
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

            
            if gen.include_current_collectors
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
                % standard length (could be moved to a paramobj geometrical field)
                xlength = 1e-6*[25, 64, 15, 57, 15]';
            else
                nxs = [nenx; sepnx; penx];
                xlength = 1e-6*[64, 15, 57]';
            end
            
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
            if gen.include_current_collectors
                gen.ccnenx = gen.ccnenx*fac;
                gen.ccpenx = gen.ccpenx*fac;
            end
            
        end
            
        function paramobj = setupElectrolyte(gen, paramobj, params)
            
            if gen.include_current_collectors
                params.cellind = gen.ccnenx + (1 : (gen.nenx + gen.sepnx + gen.penx))';
                params.Separator.cellind = gen.ccnenx + gen.nenx + (1 : gen.sepnx)';
            else
                params.cellind = (1 : (gen.nenx + gen.sepnx + gen.penx))';
                params.Separator.cellind = gen.nenx + (1 : gen.sepnx)';
            end
            paramobj = setupElectrolyte@BatteryGenerator(gen, paramobj, params);
        end
        
        function paramobj = setupElectrodes(gen, paramobj, params)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            cc  = 'CurrentCollector';
            am = 'ActiveMaterial';
            
            sepnx  = gen.sepnx; 
            nenx   = gen.nenx;
            penx   = gen.penx;
            
            if gen.include_current_collectors
                
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
            
            else
                
                %% parameters for negative electrode

                params.(ne).cellind = (1 : nenx)';
                params.(ne).(am).cellind = (1 : nenx)';
                
                % boundary setup for negative current collector
                params.(ne).(am).bcfaces = 1;
                params.(ne).(am).bccells = 1;
                
                %% parameters for positive electrode
                
                pe_indstart = nenx + sepnx;
                params.(pe).cellind =  pe_indstart + (1 :  penx)';
                params.(pe).(am).cellind = pe_indstart + (1 : penx)';
                
                % boundary setup for positive current collector
                params.(pe).(am).bcfaces = penx + 1;
                params.(pe).(am).bccells = penx;
                
            end
            
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
