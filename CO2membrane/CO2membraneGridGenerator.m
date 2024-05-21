classdef CO2membraneGridGenerator
    
    properties


        Feed     % structure with following fields
                 % - facearea
        Permeate % structure with following fields
                 % - facearea
        
        nx       % Discretization number
        length   % Length
        
    end
    
    methods

        function gen = CO2membraneGridGenerator()

            gen = gen.setupDefault();

        end

        function gen = setupDefault(gen)
        % Setup some default values

            gen.nx                = 10;
            gen.length            = 1;
            gen.Feed.facearea     = 1;
            gen.Permeate.facearea = 1;
            
        end
        
        function [inputparams, gen] = updateInputParams(gen, inputparams)

            [inputparams, gen] = gen.setupGrids(inputparams);
            [inputparams, gen] = gen.setupExternalCouplingTerms(inputparams);
            
        end

        function [inputparams, gen] = setupGrids(gen, inputparams)

            components = {'Feed', 'Permeate'};

            for icomp = 1 : numel(components)
                comp = components{icomp};
                clear params
                params.facearea = gen.(comp).facearea;
                [inputparams.(comp), gen] = setupGrid(gen, inputparams.(comp), params)
            end
            
        end
        
        function [inputparams, gen] = setupExternalCouplingTerms(gen, inputparams)

            components = {'Feed', 'Permeate'};

            for icomp = 1 : numel(components)
                comp = components{icomp};
                [inputparams.(comp), gen] = setupExternalCouplingTerm(gen, inputparams.(comp), params)
            end
            
        end
        
        function [inputparams, gen] = setupGrid(gen, inputparams, params)

            % Cartesian Grid 
            G = cartGrid(gen.nx, gen.length);
            
            % Setup parent grid with given face area
            parentGrid = Grid(G, 'faceArea', params.faceArea);

            % The component subgrid is the whole grid in this case
            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');

            inputparams.G = G;
            
        end
        
        function [inputparams, gen] = setupExternalCouplingTerm(gen, inputparams)

            bcfaces = [1;
                       gen.nx + 1];
            
            bccells = [1;
                       gen.nx];

            compnames = {'inner'};
            coupTerm = couplingTerm('boundary coupling', compnames);
            coupTerm.couplingfaces = bcfaces;
            coupTerm.couplingcells = bccells;

        end
        
    end
    
end

