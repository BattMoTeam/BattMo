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

            error('Virtual Method')
            
        end


    end
    
    
end
