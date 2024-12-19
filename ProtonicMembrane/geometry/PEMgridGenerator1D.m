classdef PEMgridGenerator1D < PEMgridGenerator


    properties

        % Length
        xlength
        % Discretization number
        N
        % Section area
        faceArea 
        
    end
    
    methods

        function gen = PEMgridGenerator1D()
            
            gen = gen@PEMgridGenerator();
            
        end

        function [inputparams, gen] = updateInputParams(gen, inputparams, params)

            inputparams = gen.setupInputParams(inputparams, []);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
            
            G = cartGrid(gen.N, gen.xlength);

            parentGrid = Grid(G, 'faceArea', gen.faceArea);
            
            G = genSubGrid(parentGrid, (1 : parentGrid.getNumberOfCells())');
            
            inputparams.G = G;
            gen.parentGrid = parentGrid;
            
        end

        function inputparams = setupElectrodeElectrolyteCoupTerm(gen, inputparams, params)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';

            G = gen.parentGrid;
            
            params.(an).couplingcells = 1;
            params.(an).couplingfaces = 1;
            params.(ct).couplingcells = G.getNumberOfCells();
            params.(ct).couplingfaces = G.getNumberOfFaces();

            inputparams = setupElectrodeElectrolyteCoupTerm@PEMgridGenerator(gen, inputparams, params);
            
        end

        function inputparams = setupControlArea(gen, inputparams)
        % Here inputparams is instance of ProtonicMembraneControlInputParams
            
            inputparams.area = gen.faceArea;
            
        end
        
    end
    
    
end
