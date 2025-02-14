classdef BatteryGenerator1DSwelling < BatteryGenerator
% setup 1D grid 
    properties

        %
        % vector of component lengths
        %
        % - x(1) : length of negative current collector (default = 25 micro meter)
        % - x(2) : length of negative active material (default = 64 micro meter)
        % - x(3) : length of separator (default = 15 micro meter)
        % - x(4) : length of positive active material (default = 57 micro meter)
        % - x(5) : length of positive current collector (default = 15 micro meter)
        % 
        xlength = 1e-6*[25; 25; 25; 71; 15];

        sepnx  = 10; % discretization number for negative current collector (default = 10)
        nenx   = 10; % discretization number for negative active material (default = 10)
        penx   = 10; % discretization number for separator (default = 10)
        nenr   = 10; % discretization for solid diffusion (default = 10)
        penr   = 10; % discretization for solid diffusion (default = 10)
        ccnenx = 10; % discretization number for positive current collector (default = 10)
        ccpenx = 10; % discretization number for positive active material (default = 10)

        %
        % refinement factor (can be used to easily increase discretization refinement)
        % see applyResolutionFactors method
        %
        fac = 1;

        % boolean : true if grid for current collectors should be included
        include_current_collectors
        % boolean : true if grid for thermal model should be included
        use_thermal

        % Face area in the transversal direction (default = 1)
        faceArea = 1; 
        
    end
    
    methods
        
        function gen = BatteryGenerator1DSwelling()
            gen = gen@BatteryGenerator();  
        end
            
        function [inputparams, gen] = updateBatteryInputParams(gen, inputparams)

            gen.include_current_collectors = inputparams.include_current_collectors;
            gen.use_thermal = inputparams.use_thermal;
            inputparams = gen.setupBatteryInputParams(inputparams, []);

            % We define some shorthand names for simplicity.
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            sep     = 'Separator';
            thermal = 'ThermalModel';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';

            % we update all the grids to adjust grid to faceArea
            inputparams.G               = gen.adjustGridToFaceArea(inputparams.G);
            inputparams.(elyte).G       = gen.adjustGridToFaceArea(inputparams.(elyte).G);
            inputparams.(elyte).(sep).G = gen.adjustGridToFaceArea(inputparams.(elyte).(sep).G);
            
            eldes = {ne, pe};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                inputparams.(elde).(am).G = gen.adjustGridToFaceArea(inputparams.(elde).(am).G);
                if gen.include_current_collectors
                    inputparams.(elde).(cc).G = gen.adjustGridToFaceArea(inputparams.(elde).(cc).G);
                end
            end
            
            if gen.use_thermal
                inputparams.(thermal).G = gen.adjustGridToFaceArea(inputparams.(thermal).G);                
            end
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, ~)
            
            sepnx  = gen.sepnx;
            nenx   = gen.nenx;
            penx   = gen.penx;
            ccnenx = gen.ccnenx;
            ccpenx = gen.ccpenx;

            xlength = gen.xlength;

            if gen.include_current_collectors
                nxs = [ccnenx; nenx; sepnx; penx; ccpenx];
            else
                nxs = [nenx; sepnx; penx];
                xlength = xlength(2 : 4);
            end
            
            x = xlength./nxs;
            x = rldecode(x, nxs);
            x = [0; cumsum(x)];

            G = tensorGrid(x);
            G = computeGeometry(G); 
            
            parentGrid = Grid(G, 'faceArea', gen.faceArea);

            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G  = G;
            gen.parentGrid = parentGrid;
            
        end

        function gen = applyResolutionFactors(gen)
            
            fac = gen.fac;
            
            gen.sepnx  = gen.sepnx*fac;
            gen.nenx   = gen.nenx*fac;
            gen.penx   = gen.penx*fac;
            gen.ccnenx = gen.ccnenx*fac;
            gen.ccpenx = gen.ccpenx*fac;
            
        end
            
        function inputparams = setupElectrolyte(gen, inputparams, params)
            
            if gen.include_current_collectors
                params.cellind = gen.ccnenx + (1 : (gen.nenx + gen.sepnx + gen.penx))';
                params.Separator.cellind = gen.ccnenx + gen.nenx + (1 : gen.sepnx)';
            else
                params.cellind = (1 : (gen.nenx + gen.sepnx + gen.penx))';
                params.Separator.cellind = gen.nenx + (1 : gen.sepnx)';
            end
            inputparams = setupElectrolyte@BatteryGenerator(gen, inputparams, params);
        end
        
        function inputparams = setupElectrodes(gen, inputparams, params)

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
            
            inputparams = setupElectrodes@BatteryGenerator(gen, inputparams, params);

        end
                
        function inputparams = setupThermalModel(gen, inputparams, params)

            params.couplingfaces = [];
            params.couplingcells = (1 : gen.G.cells.num)';
            inputparams = setupThermalModel@BatteryGenerator(gen, inputparams, params);
            
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
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
