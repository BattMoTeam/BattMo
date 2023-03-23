classdef ElectrolyserGridGenerator


    properties
        % Global grid
        % It is stored here because it is shared by many setup functions
        G

    end

    methods

        function [paramobj, gen] = updateElectrolyserInputParams(gen, paramobj, params)
        % this function is the main class function as it returns an updated paramobj object with grid structure
            error('virtual function');
        end

        function [paramobj, gen] = setupElectrolyserInputParams(gen, paramobj, params)
        % main function : add grid and coupling to paramobj structure
            inm = 'IonomerMembrane';

            [paramobj, gen] = gen.setupGrid(paramobj, params);
            paramobj.(inm)  = gen.setupIonomerGrid(paramobj.(inm), params.(inm));
            paramobj = gen.setupEvolutionElectrodes(paramobj, params);
            paramobj = setupElectrodeIonomerCoupTerm(gen, paramobj, params);

        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
        % paramobj is instance of BatteryInputParams
            error('virtual function');
        end

        function paramobj = setupIonomerGrid(gen, paramobj, params)

            paramobj.G = genSubGrid(gen.G, params.cellind);

        end

        function paramobj = setupEvolutionElectrodes(gen, paramobj, params)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';

            % We add the electrode type to the params structure (This information can be used by derived classes and may
            % then simplify setup)
            params.(oer).electrodeType = 'Oxygen';
            params.(her).electrodeType = 'Hydrogen';

            paramobj.(oer) = gen.setupEvolutionElectrode(paramobj.(oer), params.(oer));
            paramobj.(her) = gen.setupEvolutionElectrode(paramobj.(her), params.(her));

        end

        function paramobj = setupEvolutionElectrode(gen, paramobj, params)

            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';

            paramobj.G       = genSubGrid(gen.G, params.cellind);
            paramobj.(ptl).G = genSubGrid(gen.G, params.(ptl).cellind);
            paramobj.(ctl).G = genSubGrid(gen.G, params.(ctl).cellind);

            if paramobj.(ctl).include_dissolution
                dm = 'DissolutionModel';
                paramobj.(ctl).(dm).G = paramobj.(ctl).G;
            end

            %  setup coupling between catalyst layer and porous transport layer;
            G_ptl = paramobj.(ptl).G;
            G_ctl = paramobj.(ctl).G;
            G     = G_ctl.mappings.parentGrid;
            cells2 = (1 : G_ctl.cells.num)'; % all the cells of the catalyst layer are coupled to the porous transport layer
            
            mapping = zeros(G.cells.num, 1);
            mapping(G_ptl.mappings.cellmap) = (1 : G_ptl.cells.num)'; % mapping from parent grid to ptl grid;

            pcells = G_ctl.mappings.cellmap(cells2); % ctl cells indexed in parent grid
            cells1 = mapping(pcells); % ptl cells that are coupled.

            compnames = {ptl, ctl};
            coupTerm = couplingTerm(sprintf('%s-%s', ptl, ctl), compnames);
            coupTerm.couplingcells = [cells1, cells2];

            paramobj.couplingTerm = coupTerm;
            
            % setup external coupling

            compnames = {ptl};
            coupTerm = couplingTerm(sprintf('bc-%s', ptl), compnames);

            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;

            paramobj.(ptl).externalCouplingTerm = coupTerm;

        end


        function paramobj = setupElectrodeIonomerCoupTerm(gen, paramobj, params)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ctl = 'CatalystLayer';
            inm = 'IonomerMembrane';

            couplingTerms = {};

            G     = paramobj.G;
            G_inm = paramobj.(inm).G;

            mapping = zeros(G.cells.num, 1);
            mapping(G_inm.mappings.cellmap) = (1 : G_inm.cells.num)';

            eldes = {oer, her};

            % We setup also ionomer tortuosity
            tortuosity = ones(G_inm.cells.num, 1);
            
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                G_ctl = paramobj.(elde).(ctl).G;

                G = G_ctl.mappings.parentGrid;

                cells1 = (1 : G_ctl.cells.num)';
                pcells = G_ctl.mappings.cellmap(cells1);

                cells2 = mapping(pcells);

                compnames = {elde, inm};
                coupTerm = couplingTerm(sprintf('%s-%s', elde, inm), compnames);
                coupTerm.couplingcells =  [cells1, cells2];
                coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty

                couplingTerms{end + 1} = coupTerm;

                tortuosity(cells2) = paramobj.(elde).(ctl).tortuosity;
                
            end

            paramobj.couplingTerms = couplingTerms;
            paramobj.(inm).tortuosity = tortuosity;
            

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
