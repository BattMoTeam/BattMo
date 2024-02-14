classdef GasSupplyGridGenerator

    properties
        
        parentGrid

    end


    methods
        
        function [inputparams, gen] = updateInputParams(gen, inputparams, params)
            
            error('virtual function');
            
        end
        
        function [inputparams, gen] = setupInputParams(gen, inputparams, params)
        % This is the main function which should be called by the derived class at the end of the :meth:`updateInputParams` method

            
            [inputparams, gen] = gen.setupGrid(inputparams, params);
            [inputparams, gen] = gen.setupExternalCoupling(inputparams, params);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
            
            error('virtual function');
            
        end

        function [inputparams, gen] = setupExternalCoupling(gen, inputparams, params)
            
            error('virtual function');
            
        end


    end


end
