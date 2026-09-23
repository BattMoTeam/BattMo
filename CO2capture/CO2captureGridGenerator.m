classdef CO2captureGridGenerator
    
    properties

        nx       % Discretization number
        length   % Length
        areas
        
    end
    
    methods

        function gen = CO2captureGridGenerator()

            gen = gen.setupDefault();

        end

        function gen = setupDefault(gen)
        % Setup some default values

            gen.nx     = 100;
            gen.length = 1;

            gen.areas = struct('feed', 1, ...
                               'permeate', 1);
            
        end
        
        function [inputparams, gen] = updateInputParams(gen, inputparams)

            [inputparams, gen] = gen.setupGrids(inputparams);
            [inputparams, gen] = gen.setupCrossCouplingTerm(inputparams);
            [inputparams, gen] = gen.setupControlCouplingTerms(inputparams);
            
        end

        function [inputparams, gen] = setupGrids(gen, inputparams)

            components = {'Feed', 'Permeate'};

            for icomp = 1 : numel(components)
                comp = components{icomp};
                switch comp
                  case 'Feed'
                    area = gen.areas.feed;
                  case 'Permeate'
                    area = gen.areas.permeate;
                  otherwise
                    error('Unknown component: %s', comp);
                end
                [inputparams.(comp), gen] = setupGrid(gen, inputparams.(comp), area);
            end
            
        end
        
        function [inputparams, gen] = setupControlCouplingTerms(gen, inputparams)

            components = {'Feed', 'Permeate'};

            for icomp = 1 : numel(components)
                comp = components{icomp};
                [inputparams.(comp), gen] = setupControlCouplingTerm(gen, inputparams.(comp));
            end
            
        end
        
        
        function [inputparams, gen] = setupGrid(gen, inputparams, area)

            % Cartesian Grid 
            G = cartGrid(gen.nx, gen.length);
            
            % Setup parent grid with given face area
            parentGrid = Grid(G, 'faceArea', area);

            % The component subgrid is the whole grid in this case
            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G = G;
            
        end
        
        function [inputparams, gen] = setupCrossCouplingTerm(gen, inputparams)

            
            bccells = [(1 : gen.nx)', (1 : gen.nx)'];

            compnames = {'Feed', 'Permeate'};
            coupTerm = couplingTerm('Feed-Permeate', compnames);

            coupTerm.couplingcells = bccells;

            inputparams.couplingTerms{end + 1} = coupTerm;
            
        end

        
        function [inputparams, gen] = setupControlCouplingTerm(gen, inputparams, comp)

            bcfaces = [1; gen.nx + 1];
            bccells = [1; gen.nx];

            coupname = 'boundary faces';
            compnames = {'channel'};
            coupTerm = couplingTerm(coupname, compnames);
            coupTerm.couplingfaces = bcfaces;
            coupTerm.couplingcells = bccells;

            inputparams.couplingTerms{end + 1} = coupTerm;

            coupname = 'control';
            compnames = {'channel boundary faces', 'control'};
            coupTerm = couplingTerm(coupname, compnames);
            coupTerm.couplingfaces = [bcfaces, [1; 2]];
            
            inputparams.couplingTerms{end + 1} = coupTerm;

        end

    end
    
end

