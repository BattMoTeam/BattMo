classdef ElectrolyserGridGenerator


    properties
        % Global grid
        % It is stored here because it is shared by many setup functions
        G

    end

    methods

        function [inputparams, gen] = updateElectrolyserInputParams(gen, inputparams, params)
        % this function is the main class function as it returns an updated inputparams object with grid structure
            error('virtual function');
        end

        function [inputparams, gen] = setupElectrolyserInputParams(gen, inputparams, params)
        % main function : add grid and coupling to inputparams structure
            inm = 'IonomerMembrane';

            [inputparams, gen] = gen.setupGrid(inputparams, params);
            inputparams.(inm)  = gen.setupIonomerGrid(inputparams.(inm), params.(inm));
            inputparams = gen.setupEvolutionElectrodes(inputparams, params);
            inputparams = setupElectrodeIonomerCoupTerm(gen, inputparams, params);

        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
        % inputparams is instance of BatteryInputParams
            error('virtual function');
        end

        function inputparams = setupIonomerGrid(gen, inputparams, params)

            inputparams.G = genSubGrid(gen.G, params.cellind);

        end

        function inputparams = setupEvolutionElectrodes(gen, inputparams, params)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';

            % We add the electrode type to the params structure (This information can be used by derived classes and may
            % then simplify setup)
            params.(oer).electrode_type = 'Oxygen';
            params.(her).electrode_type = 'Hydrogen';

            inputparams.(oer) = gen.setupEvolutionElectrode(inputparams.(oer), params.(oer));
            inputparams.(her) = gen.setupEvolutionElectrode(inputparams.(her), params.(her));

        end

        function inputparams = setupEvolutionElectrode(gen, inputparams, params)

            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';

            inputparams.G       = genSubGrid(gen.G, params.cellind);
            inputparams.(ptl).G = genSubGrid(gen.G, params.(ptl).cellind);
            inputparams.(ctl).G = genSubGrid(gen.G, params.(ctl).cellind);

            if inputparams.(ctl).include_dissolution
                dm = 'DissolutionModel';
                inputparams.(ctl).(dm).G = inputparams.(ctl).G;
            end

            %  setup coupling between catalyst layer and porous transport layer;
            G_ptl = inputparams.(ptl).G;
            G_ctl = inputparams.(ctl).G;
            G     = G_ctl.mappings.parentGrid;
            cells2 = (1 : G_ctl.cells.num)'; % all the cells of the catalyst layer are coupled to the porous transport layer

            mapping = zeros(G.getNumberOfCells(), 1);
            mapping(G_ptl.mappings.cellmap) = (1 : G_ptl.cells.num)'; % mapping from parent grid to ptl grid;

            pcells = G_ctl.mappings.cellmap(cells2); % ctl cells indexed in parent grid
            cells1 = mapping(pcells); % ptl cells that are coupled.

            compnames = {ptl, ctl};
            coupTerm = couplingTerm(sprintf('%s-%s', ptl, ctl), compnames);
            coupTerm.couplingcells = [cells1, cells2];

            inputparams.couplingTerm = coupTerm;

            % setup external coupling

            compnames = {ptl};
            coupTerm = couplingTerm(sprintf('bc-%s', ptl), compnames);

            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;

            inputparams.(ptl).externalCouplingTerm = coupTerm;

        end


        function inputparams = setupElectrodeIonomerCoupTerm(gen, inputparams, params)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ctl = 'CatalystLayer';
            inm = 'IonomerMembrane';

            couplingTerms = {};

            G     = inputparams.G;
            G_inm = inputparams.(inm).G;

            mapping = zeros(G.getNumberOfCells(), 1);
            mapping(G_inm.mappings.cellmap) = (1 : G_inm.cells.num)';

            eldes = {oer, her};

            % We setup also ionomer tortuosity
            tortuosity = ones(G_inm.cells.num, 1);

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                G_ctl = inputparams.(elde).(ctl).G;

                G = G_ctl.mappings.parentGrid;

                cells1 = (1 : G_ctl.cells.num)';
                pcells = G_ctl.mappings.cellmap(cells1);

                cells2 = mapping(pcells);

                compnames = {elde, inm};
                coupTerm = couplingTerm(sprintf('%s-%s', elde, inm), compnames);
                coupTerm.couplingcells =  [cells1, cells2];
                coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty

                couplingTerms{end + 1} = coupTerm;

                tortuosity(cells2) = inputparams.(elde).(ctl).tortuosity;

            end

            inputparams.couplingTerms = couplingTerms;
            inputparams.(inm).tortuosity = tortuosity;


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
