classdef PEMgridGenerator1D < PEMgridGenerator


    properties

        xlength  = 22*micro*meter;
        N        = 100;
        faceArea = 1;
        
    end
    
    methods

        function gen = PEMgridGenerator1D()
            
            gen = gen@PEMgridGenerator();
            
        end

        function [paramobj, gen] = updatePEMinputParams(gen, paramobj, params)

            paramobj = gen.setupPEMinputParams(paramobj, []);
            
            elyte = 'Electrolyte';
            
            paramobj.G         = gen.adjustGridToFaceArea(paramobj.G);
            paramobj.(elyte).G = gen.adjustGridToFaceArea(paramobj.(elyte).G);
            
            paramobj.dx    = gen.xlength/gen.N;
            paramobj.faceArea = gen.faceArea;
            
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
            
            G = cartGrid(gen.N, gen.xlength);
            G = computeGeometry(G);
            
            paramobj.G = G;
            gen.G = G;
            
        end
        
        function paramobj = setupElectrolyte(gen, paramobj, params)
        % Method that setups the grid and the coupling for the electrolyte model

            params.cellind = (1 : gen.N)';
            paramobj = setupElectrolyte@PEMgridGenerator(gen, paramobj, params);
            
        end

        function paramobj = setupElectrodeElectrolyteCoupTerm(gen, paramobj, params)

            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            G = paramobj.Electrolyte.G;
            
            couplingTerms = {};
            
            coupterm = couplingTerm('Anode-Electrolyte', {an, elyte});
            coupterm.couplingcells = [1, 1];
            coupterm.couplingfaces = [1, 1];
            couplingTerms{end + 1} = coupterm;
            
            coupterm = couplingTerm('Cathode-Electrolyte', {an, elyte});
            coupterm.couplingcells = [1, G.cells.num];
            coupterm.couplingfaces = [1, G.faces.num];
            couplingTerms{end + 1} = coupterm;

            paramobj.couplingTerms = couplingTerms;
            
        end

        function G = adjustGridToFaceArea(gen, G);

            fa = gen.faceArea;

            G.faces.areas   = fa*G.faces.areas;
            G.faces.normals = fa*G.faces.normals;
            G.cells.volumes = fa*G.cells.volumes;

        end
        
    end
    
    
end
