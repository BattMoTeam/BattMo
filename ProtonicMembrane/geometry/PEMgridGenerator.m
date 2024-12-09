classdef PEMgridGenerator

    properties

        parentGrid
        
    end

    methods
        
        function [inputparams, gen] = updateInputParams(gen, inputparams, params)
            
            error('virtual function');
            
        end
        
        function [inputparams, gen] = setupInputParams(gen, inputparams, params)
        % This is the main function which should be called by the derived class at the end of the :meth:`updateInputParams` method
        % The function set up the grid and the coupling terms and add those in the :code:`inputparams` structure which is
        % an instance of :class:`InputParams <PEM.PEMinputParams>`
        % 
            
            [inputparams, gen] = gen.setupGrid(inputparams, params);

            params.Electrolyte.cellinds = inputparams.G.mappings.cellmap;
            inputparams.Electrolyte = gen.setupElectrolyte(inputparams.Electrolyte, params.Electrolyte);

            inputparams = gen.setupElectrodeElectrolyteCoupTerm(inputparams);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
            
            error('virtual function');
            
        end
        
        function inputparams = setupElectrolyte(gen, inputparams, params)
        % Method that setups the grid and the coupling for the electrolyte model
        % Here inputparams is instance of ProtonicMembraneElectrolyteInputParams
           
            inputparams = gen.setupElectrolyteGrid(inputparams, params);
            
        end

        function inputparams = setupElectrolyteGrid(gen, inputparams, params)
        % Setup the grid for the electrolyte
            
            % Default setup
            inputparams.G = genSubGrid(gen.parentGrid, params.cellinds);
            
        end

        function inputparams = setupElectrodeElectrolyteCoupTerm(gen, inputparams, params)


            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            params.(an).couplingcells;
            params.(ct).couplingcells;
            
            couplingTerms = {};
            
            coupterm = couplingTerm('Anode-Electrolyte', {an, elyte});
            nc = numel(params.(an).couplingcells);
            coupterm.couplingcells = [(1 : nc)', params.(an).couplingcells];
            coupterm.couplingfaces = [(1 : nc)', params.(an).couplingfaces];
            couplingTerms{end + 1} = coupterm;
            
            coupterm = couplingTerm('Cathode-Electrolyte', {ct, elyte});
            nc = numel(params.(ct).couplingcells);
            coupterm.couplingcells = [(1 : nc)', params.(ct).couplingcells];
            coupterm.couplingfaces = [(1 : nc)', params.(ct).couplingfaces];
            couplingTerms{end + 1} = coupterm;

            inputparams.couplingTerms = couplingTerms;
            
        end


    end
    
    
end
