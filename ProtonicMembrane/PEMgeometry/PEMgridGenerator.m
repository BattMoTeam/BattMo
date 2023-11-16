classdef PEMgridGenerator

    properties

        G
        
    end

    methods
        
        function [paramobj, gen] = updatePEMinputParams(gen, paramobj, params)
            
            error('virtual function');
            
        end
        
        function [paramobj, gen] = setupPEMinputParams(gen, paramobj, params)
        % This is the main function which should be called by the derived class at the end of the :meth:`updatePEMinputParams` method
        % The function set up the grid and the coupling terms and add those in the :code:`paramobj` structure which is
        % an instance of :class:`PEMinputParams <PEM.PEMinputParams>`
        % 
            
            [paramobj, gen] = gen.setupGrid(paramobj, params);

            params = pickField(params, 'Electrolyte');
            paramobj.Electrolyte = gen.setupElectrolyte(paramobj.Electrolyte, params);

            paramobj = gen.setupElectrodeElectrolyteCoupTerm(paramobj);
            
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
            
            error('virtual function');
            
        end
        
        function paramobj = setupElectrolyte(gen, paramobj, params)
        % Method that setups the grid and the coupling for the electrolyte model
        % Here paramobj is instance of ProtonicMembraneElectrolyteInputParams
           
            paramobj = gen.setupElectrolyteGrid(paramobj, params);
            
        end

        function paramobj = setupElectrolyteGrid(gen, paramobj, params)
        % Setup the grid for the electrolyte
            
            % Default setup
            paramobj.G = genSubGrid(gen.G, params.cellind);
            
        end

        function paramobj = setupElectrodeElectrolyteCoupTerm(gen, paramobj, params)


            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            G = paramobj.Electrolyte.G;

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

            paramobj.couplingTerms = couplingTerms;
            
        end


    end
    
    
end
