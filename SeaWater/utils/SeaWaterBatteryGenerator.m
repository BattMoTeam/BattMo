classdef SeaWaterBatteryGenerator < BatteryGenerator

    methods

        function gen = SeaWaterBatteryGenerator()
            gen = gen@BatteryGenerator();
        end

        function inputparams = setupBatteryInputParams(gen, inputparams, params)
            [inputparams, gen] = gen.setupGrid(inputparams, params);
            inputparams.Electrolyte = gen.setupElectrolyte(inputparams.Electrolyte, params);
            inputparams = gen.setupElectrodes(inputparams, params);
            inputparams = gen.setupElectrodeElectrolyteCoupTerm(inputparams);
        end

        function inputparams = setupElectrolyte(gen, inputparams, params)
            inputparams = gen.setupElectrolyteGrid(inputparams, params);
        end


        function inputparams = setupElectrodes(gen, inputparams, params)

            ctde = 'Cathode';
            ande = 'Anode';
            % setup Cathode
            inputparams.(ctde) = gen.setupElectrodeGrid(inputparams.(ctde), params.(ctde));
            inputparams.(ctde) = gen.setupElectrodeBcCoupTerm(inputparams.(ctde), params.(ctde));
            % setup Anode
            inputparams.(ande) = gen.setupElectrodeGrid(inputparams.(ande), params.(ande));
            inputparams.(ande) = gen.setupElectrodeBcCoupTerm(inputparams.(ande), params.(ande));

        end

        function inputparams = setupElectrodeBcCoupTerm(gen, inputparams, params)

            % default setup
            compname = 'Electrode';
            compnames = {compname};
            coupname = sprintf('Exterior-%s', compname);
            coupTerm = couplingTerm(coupname, compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;

            inputparams.externalCouplingTerm = coupTerm;

        end

        function inputparams = setupElectrodeElectrolyteCoupTerm(gen, inputparams, params)
        % inputparams is instance of SeaWaterBatteryInputParams
        % setup inputparams.couplingTerms
            ctde  = 'Cathode';
            ande  = 'Anode';
            elyte = 'Electrolyte';

            couplingTerms = {};

            G_ctde = inputparams.(ctde).G;
            G_elyte = inputparams.(elyte).G;

            % parent Grid
            G = G_ctde.parentGrid;

            % All the cells from Cathode are coupled with Electrolyte
            cells1 = (1 : G_ctde.getNumberOfCells())';
            pcells = G_ctde.mappings.cellmap(cells1);

            mapping = zeros(G.getNumberOfCells(), 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.getNumberOfCells())';
            cells2 = mapping(pcells);

            compnames = {'Cathode', 'Electrolyte'};
            coupTerm = couplingTerm('Cathode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty

            couplingTerms{end + 1} = coupTerm;

            G_ande = inputparams.(ande).G;
            G_elyte = inputparams.(elyte).G;

            % parent Grid
            G = G_ande.parentGrid;

            % All the cells from Anode are coupled with Electrolyte
            cells1 = (1 : G_ande.getNumberOfCells())';
            pcells = G_ande.mappings.cellmap(cells1);

            mapping = zeros(G.getNumberOfCells(), 1);
            mapping(G_elyte.mappings.cellmap) = (1 : G_elyte.getNumberOfCells())';
            cells2 = mapping(pcells);

            compnames = {'Anode', 'Electrolyte'};
            coupTerm = couplingTerm('Anode-Electrolyte', compnames);
            coupTerm.couplingcells =  [cells1, cells2];
            coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty

            couplingTerms{end + 1} = coupTerm;

            inputparams.couplingTerms = couplingTerms;

        end

        function inputparams = setupCurrentCollectorElectrodeActiveComponentCoupTerm(gen, inputparams, params)
        % inputparams is instance of ElectrodeInputParams
        % setup inputparams.couplingTerm
            eac = 'ElectrodeActiveComponent';
            cc  = 'CurrentCollector';

            compnames = {'CurrentCollector', 'ElectrodeActiveComponent'};
            coupTerm = couplingTerm('CurrentCollector-ElectrodeActiveComponent', compnames);

            G_eac = inputparams.(eac).G;
            G_cc = inputparams.(cc).G;

            cctbl.faces = (1 : G_cc.getNumberOfFaces())';
            cctbl.globfaces = G_cc.mappings.facemap;
            cctbl = IndexArray(cctbl);

            eactbl.faces = (1 : G_eac.getNumberOfFaces())';
            eactbl.globfaces = G_eac.mappings.facemap;
            eactbl = IndexArray(eactbl);

            gen = CrossIndexArrayGenerator();
            gen.tbl1 = cctbl;
            gen.tbl2 = eactbl;
            gen.replacefds1 = {{'faces', 'faces1'}};
            gen.replacefds2 = {{'faces', 'faces2'}};
            gen.mergefds = {'globfaces'};
            tbl = gen.eval();

            cc_coupfaces = tbl.get('faces1');
            eac_coupfaces = tbl.get('faces2');

            cc_coupcells = sum(G_cc.faces.neighbors(cc_coupfaces, :), 2);
            eac_coupcells = sum(G_eac.faces.neighbors(eac_coupfaces, :), 2);

            coupTerm.couplingfaces = [cc_coupfaces, eac_coupfaces];
            coupTerm.couplingcells = [cc_coupcells, eac_coupcells];

            inputparams.couplingTerm = coupTerm;

        end

        function inputparams = setupCurrentCollectorBcCoupTerm(gen, inputparams, params)
        % inputparams is instance of CurrentCollectorInputParams
        % setup inputparams.couplingTerm

            % default setup
            compnames = {'CurrentCollector'};
            coupTerm = couplingTerm('bc-CurrentCollector', compnames);
            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;

            inputparams.couplingTerm = coupTerm;

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
