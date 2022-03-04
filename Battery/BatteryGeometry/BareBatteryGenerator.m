classdef BareBatteryGenerator
% Object that add grids and coupling terms to a paramobj instance of BareBatteryInputParams (through method 
% updateBatteryInputParams)
%
% This class goes through the whole setup and is meant to be used as a base class.
%
% This class (or rather an subclass of it) can also be used to setup values of the paramobj instance that is sent to it
% (just overload the updateBatteryInputParams method by adding the desired setup)
%
% Example : BareBatteryGenerator1D.m

    properties
        % Global grid 
        % It is stored here because is shared by many setup functions
        % It is setup by function setupBatteryGrid
        G
    end

    methods
        
        function paramobj = updateBatteryInputParams(gen, paramobj)
            error('virtual function - should call setupBatteryInputParams with some argument for params');
        end
        
        function paramobj = setupBatteryInputParams(gen, paramobj, params)
        % main function : add grid and coupling to paramobj structure
            [paramobj, gen] = gen.setupGrid(paramobj, params);
            paramobj.Electrolyte = gen.setupElectrolyte(paramobj.Electrolyte, params);
            paramobj = gen.setupElectrodes(paramobj, params);
            paramobj = gen.setupElectrodeElectrolyteCoupTerm(paramobj);
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
        % setup paramobj.G and update gen.G 
            error('virtual function');
        end
        
        function paramobj = setupElectrolyte(gen, paramobj, params)
        % paramobj is instance of ElectrolyteInputParams
            paramobj = gen.setupElectrolyteGrid(paramobj, params);
            paramobj.Separator = gen.setupSeparatorGrid(paramobj.Separator, params.Separator);
        end

        function paramobj = setupElectrolyteGrid(gen, paramobj, params)
        % paramobj is instance of ElectrolyteInputParams
        % setup paramobj.G

            % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
            
        end
        
        function paramobj = setupSeparatorGrid(gen, paramobj, params)
        % paramobj is instance of SeparatorInputParams
        % setup paramobj.G
            
           % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
            
        end

        function paramobj = setupElectrodes(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            % setup Negative Electrode 
            paramobj.(ne) = gen.setupElectrodeGrid(paramobj.(ne), params.(ne));
            % setup Negative Electrode Exterior coupling
            params.(ne).compname = ne;
            paramobj = gen.setupElectrodeBcCoupTerm(paramobj, params.(ne));
            % setup Positive Electrode
            paramobj.(pe) = gen.setupElectrodeGrid(paramobj.(pe), params.(pe));
            % setup Positive Electrode Exterior coupling
            params.(pe).compname = pe;
            paramobj = gen.setupElectrodeBcCoupTerm(paramobj, params.(pe));            
        end
                
        
        function paramobj = setupElectrodeGrid(gen, paramobj, params)
        % paramobj is instance of ActiveMaterialInputParams
        % setup paramobj.G
            
            % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
        end

        function paramobj = setupElectrodeElectrolyteCoupTerm(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
        % setup paramobj.couplingTerms
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            
            couplingTerms = {};
            
            G_ne = paramobj.(ne).G;
            G_elyte = paramobj.(elyte).G;
            
            % parent Grid
            G = G_ne.mappings.parentGrid;
            
            % All the cells from NegativeElectrode are coupled with Electrolyte
            cells1 = (1 : G_ne.cells.num)';
            pcells = G_ne.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'NegativeElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('NegativeElectrode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
            couplingTerms{end + 1} = coupTerm;
            
            G_pe = paramobj.(pe).G;
            G_elyte = paramobj.(elyte).G;
            
            % parent Grid
            G = G_pe.mappings.parentGrid;
            
            % All the cells from PositiveElectrode are coupled with Electrolyte
            cells1 = (1 : G_pe.cells.num)';
            pcells = G_pe.mappings.cellmap(cells1);
            
            mapping = zeros(G.cells.num, 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.cells.num)';
            cells2 = mapping(pcells);
            
            compnames = {'PositiveElectrode', 'Electrolyte'};
            coupTerm = couplingTerm('PositiveElectrode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty
            
            couplingTerms{end + 1} = coupTerm;
               
            if isempty(paramobj.couplingTerms)
                paramobj.couplingTerms = couplingTerms;
            else
                paramobj.couplingTerms = {paramobj.couplingTerms{:}, couplingTerms{:}};
            end

        end

        function paramobj = setupElectrodeBcCoupTerm(gen, paramobj, params)
        % paramobj is instance of BareBatteryInputParams
        % add element in paramobj.couplingTerms
            
            % default setup
            compname = params.compname;
            compnames = {compname};
            coupname = sprintf('Exterior-%s', compname);
            coupTerm = couplingTerm(coupname, compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;
            
            couplingTerms = paramobj.couplingTerms;
            if isempty(couplingTerms)
                paramobj.couplingTerms = {coupTerm};
            else
                paramobj.couplingTerms = {couplingTerms{:}, coupTerm};
            end
        end
        
    end
    
    
end





%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
