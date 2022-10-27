classdef BatteryGenerator1D < BatteryGenerator
% setup 1D grid 
    properties
        
        xlength = 1e-6*[25; 64; 15; 57; 15]; % length of components in x direction - default values
                                             % x(1) : length of negative current collector
                                             % x(2) : length of negative activie material
                                             % x(3) : length of separator
                                             % x(4) : length of positive current collector
                                             % x(5) : length of positive activie material

        sepnx  = 10; % discretization number for negative current collector - default value
        nenx   = 10; % discretization number for negative activie material  - default value
        penx   = 10; % discretization number for separator                  - default value
        nenr   = 10; % discretization for solid diffusion                   - default value
        penr   = 10; % discretization for solid diffusion                   - default value
        ccnenx = 10; % discretization number for positive current collector - default value
        ccpenx = 10; % discretization number for positive activie material  - default value

        % refinement factor (can be used to easily increase discretization refinement)
        % see applyResolutionFactors method
        fac = 1;

        % flag : true if grid for current collectors should be included
        include_current_collectors
        % flag : true if grid for thermal model should be included
        use_thermal

        faceArea = 1; 
        
    end
    
    methods
        
        function gen = BatteryGenerator1D()
            gen = gen@BatteryGenerator();  
        end
            
        function paramobj = updateBatteryInputParams(gen, paramobj)

            gen.include_current_collectors = paramobj.include_current_collectors;
            gen.use_thermal = paramobj.use_thermal;
            paramobj = gen.setupBatteryInputParams(paramobj, []);

            % We define some shorthand names for simplicity.
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            sep     = 'Separator';
            thermal = 'ThermalModel';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';

            % we update all the grids to adjust grid to faceArea
            paramobj.G               = gen.adjustGridToFaceArea(paramobj.G);
            paramobj.(elyte).G       = gen.adjustGridToFaceArea(paramobj.(elyte).G);
            paramobj.(elyte).(sep).G = gen.adjustGridToFaceArea(paramobj.(elyte).(sep).G);
            
            eldes = {ne, pe};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                paramobj.(elde).(am).G = gen.adjustGridToFaceArea(paramobj.(elde).(am).G);
                if gen.include_current_collectors
                    paramobj.(elde).(cc).G = gen.adjustGridToFaceArea(paramobj.(elde).(cc).G);
                end
            end
            
            if gen.use_thermal
                paramobj.(thermal).G = gen.adjustGridToFaceArea(paramobj.(thermal).G);                
            end
            
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, ~)
        % paramobj is instance of BatteryInputParams
        % setup paramobj.G
            
            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;
            ccnenx = gen.ccnenx;
            ccpenx = gen.ccpenx;

            xlength = gen.xlength;

            if gen.include_current_collectors
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
                % standard length (could be moved to a paramobj geometrical field)
            else
                nxs = [nenx; sepnx; penx];
                xlength = xlength(2 : 4);
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
            %if gen.include_current_collectors
            gen.ccnenx = gen.ccnenx*fac;
            gen.ccpenx = gen.ccpenx*fac;
            %end
            
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

        function G = adjustGridToFaceArea(gen, G);
            
            fa = gen.faceArea;

            G.faces.areas   = fa*G.faces.areas;
            G.faces.normals = fa*G.faces.normals;
            G.cells.volumes = fa*G.cells.volumes;

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
