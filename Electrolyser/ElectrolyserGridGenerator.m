classdef ElectrolyserGridGenerator


    properties
        % Global grid
        % It is stored here because is shared by many setup functions
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

            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            % We add the electrode type to the params structure (This information can be used by derived classes and may
            % then simplify setup)
            params.(her).electrodeType = 'Hydrogen';
            params.(oer).electrodeType = 'Oxygen';

            paramobj.(her) = gen.setupEvolutionElectrode(paramobj.(her), params.(her));
            paramobj.(oer) = gen.setupEvolutionElectrode(paramobj.(oer), params.(oer));

        end

        function paramobj = setupEvolutionElectrode(gen, paramobj, params)

            ptl = 'PorousTransportLayer';
            paramobj.(ptl).G = genSubGrid(gen.G, params.cellind);

            % setup external coupling

            compnames = {ptl};
            coupTerm = couplingTerm(sprintf('bc-%s', ptl), compnames);

            coupTerm.couplingfaces = params.bcfaces;
            coupTerm.couplingcells = params.bccells;

            paramobj.(ptl).externalCouplingTerm = coupTerm;

        end


        function paramobj = setupElectrodeIonomerCoupTerm(gen, paramobj, params)

            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            inm = 'IonomerMembrane';

            couplingTerms = {};

            G     = paramobj.G;
            G_inm = paramobj.(inm).G;

            mapping = zeros(G.cells.num, 1);
            mapping(G_inm.mappings.cellmap) = (1 : G_inm.cells.num)';

            eldes = {her, oer};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                G_elde = paramobj.(elde).(ptl).G;

                G = G_elde.mappings.parentGrid;

                cells1 = params.(elde).coupcellind;
                pcells = G_elde.mappings.cellmap(cells1);

                cells2 = mapping(pcells);

                compnames = {elde, inm};
                coupTerm = couplingTerm(sprintf('%s-%s', elde, inm), compnames);
                coupTerm.couplingcells =  [cells1, cells2];
                coupTerm.couplingfaces = []; % no coupling throug faces. We set it as empty

                couplingTerms{end + 1} = coupTerm;

            end

            paramobj.couplingTerms = couplingTerms;

        end


    end


end

