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

        function [inputparams, gen] = updatePEMinputParams(gen, inputparams, params)

            inputparams = gen.setupPEMinputParams(inputparams, []);
            
            elyte = 'Electrolyte';
            
            inputparams.G         = gen.adjustGridToFaceArea(inputparams.G);
            inputparams.(elyte).G = gen.adjustGridToFaceArea(inputparams.(elyte).G);
            
            inputparams.dx    = gen.xlength/gen.N;
            inputparams.faceArea = gen.faceArea;
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
            
            G = cartGrid(gen.N, gen.xlength);
            G = computeGeometry(G);
            
            inputparams.G = G;
            gen.G = G;
            
        end
        
        function inputparams = setupElectrolyte(gen, inputparams, params)
        % Method that setups the grid and the coupling for the electrolyte model

            params.cellind = (1 : gen.N)';
            inputparams = setupElectrolyte@PEMgridGenerator(gen, inputparams, params);
            
        end

        function inputparams = setupElectrodeElectrolyteCoupTerm(gen, inputparams, params)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';

            G = inputparams.G;
            
            params.(an).couplingcells = 1;
            params.(an).couplingfaces = 1;
            params.(ct).couplingcells = G.cells.num;
            params.(ct).couplingfaces = G.faces.num;

            inputparams = setupElectrodeElectrolyteCoupTerm@PEMgridGenerator(gen, inputparams, params);
            
        end

        function G = adjustGridToFaceArea(gen, G);

            fa = gen.faceArea;

            G.faces.areas   = fa*G.faces.areas;
            G.faces.normals = fa*G.faces.normals;
            G.cells.volumes = fa*G.cells.volumes;

        end
        
    end
    
    
end
