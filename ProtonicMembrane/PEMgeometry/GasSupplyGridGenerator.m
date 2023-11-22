classdef GasSupplyGridGenerator

    properties
        
        G

    end


    methods
        
        function [paramobj, gen] = updateInputParams(gen, paramobj, params)
            
            error('virtual function');
            
        end
        
        function [paramobj, gen] = setupInputParams(gen, paramobj, params)
        % This is the main function which should be called by the derived class at the end of the :meth:`updateInputParams` method

            
            [paramobj, gen] = gen.setupGrid(paramobj, params);
            [paramobj, gen] = gen.setupExternalCoupling(paramobj, params);
            
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
            
            error('virtual function');
            
        end

        function [paramobj, gen] = setupExternalCoupling(gen, paramobj, params)
            
            error('virtual function');
            
        end


    end


end
