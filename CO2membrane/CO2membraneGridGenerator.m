classdef CO2membraneGridGenerator
    
    properties

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
            
        end
        
        function [inputparams, gen] = updateInputParams(gen, inputparams)

            [inputparams, gen] = gen.setupGrids(inputparams);
            [inputparams, gen] = gen.setupExternalCouplingTerms(inputparams);
            [inputparams, gen] = gen.setupControlCouplingTerms(inputparams);
            
        end

        function [inputparams, gen] = setupGrids(gen, inputparams)

            components = {'Feed', 'Permeate'};

            for icomp = 1 : numel(components)
                comp = components{icomp};
                [inputparams.(comp), gen] = setupGrid(gen, inputparams.(comp));
            end
            
        end
        
        function [inputparams, gen] = setupExternalCouplingTerms(gen, inputparams)

            components = {'Feed', 'Permeate'};

            for icomp = 1 : numel(components)
                comp = components{icomp};
                [inputparams.(comp), gen] = setupExternalCouplingTerm(gen, inputparams.(comp));
            end
            
        end

        function [inputparams, gen] = setupControlCouplingTerms(gen, inputparams)

            components = {'Feed', 'Permeate'};

            for icomp = 1 : numel(components)
                comp = components{icomp};
                [inputparams.(comp), gen] = setupControlCouplingTerm(gen, inputparams.(comp));
            end
            
        end
        
        
        function [inputparams, gen] = setupGrid(gen, inputparams, params)

            % Cartesian Grid 
            G = cartGrid(gen.nx, gen.length);
            
            % Setup parent grid with given face area
            parentGrid = Grid(G);

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
            coupTerm = couplingTerm('boundary faces', compnames);
            coupTerm.couplingfaces = bcfaces;
            coupTerm.couplingcells = bccells;

            inputparams.couplingTerms{end + 1} = coupTerm;
            
        end

        
        function [inputparams, gen] = setupControlCouplingTerm(gen, inputparams)

            bcfaces = 1;
            bccells = 1;

            compnames = {'inner'};
            coupTerm = couplingTerm('control faces', compnames);
            coupTerm.couplingfaces = bcfaces;
            coupTerm.couplingcells = bccells;

            inputparams.couplingTerms{end + 1} = coupTerm;
            
        end

    end
    
end

