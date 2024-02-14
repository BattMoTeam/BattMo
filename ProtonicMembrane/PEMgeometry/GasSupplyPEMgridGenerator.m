classdef GasSupplyPEMgridGenerator

    properties

        parentGrid
        
        % Use to generate the structure of the internal coupling terms
        cellGridGenerator
        gasSupplyGridGenerator
        
    end

    methods
        
        function [inputparams, gen] = updateInputParams(gen, inputparams, params)
            
            error('virtual function');
            
        end
        
        function [inputparams, gen] = setupInputParams(gen, inputparams, params)
        % This is the main function which should be called by the derived class at the end of the :meth:`updateInputParams` method

            
            [inputparams, gen] = gen.setupGrid(inputparams, params);

            gen = gen.setupCellGridGenerator(inputparams, params);
            gen = gen.setupGasSupplyGridGenerator(inputparams, params);

            inputparams.Cell.G      = genSubGrid(gen.parentGrid, params.Cell.cellinds);
            inputparams.GasSupply.G = genSubGrid(gen.parentGrid, params.GasSupply.cellinds);
            
            [inputparams.Cell, gen]      = gen.setupCell(inputparams.Cell, params.Cell);
            [inputparams.GasSupply, gen] = gen.setupGasSupply(inputparams.GasSupply, params.GasSupply);

            [inputparams, gen] = gen.setupCellGasSupplyCoupling(inputparams);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
            
            error('virtual function');
            
        end
        
        function gen = setupCellGridGenerator(gen, inputparams, params)
        % setup gen.cellGridGenerator

            error('virtual function');
            
        end

        function gen = setupGasSupplyGridGenerator(gen, inputparams, params)
        % setup gen.gasSupplyGridGenerator
            
            error('virtual function');
            
        end

        function [inputparams, gen] = setupCell(gen, inputparams, params)
        % - inputparams is instance of ProtonicMembraneInputParams
        % - The setup is taken from setupInputParams method in PEMgridGenerator except for the grid setup.
        %   Note that we assume that electrolyte and PEM have same grid (this may changed...)
        % - params has the same fields as the one sent to setupInputParams method.
        %   params.Electrolyte.cellind should correspond to gen.parentGrid (which is same as gen.cellGridGenerator.parentGrid).

            cgen = gen.cellGridGenerator;
            
            inputparams.Electrolyte = cgen.setupElectrolyte(inputparams.Electrolyte, params.Electrolyte);
            
            inputparams = cgen.setupElectrodeElectrolyteCoupTerm(inputparams);
            
        end

        function [inputparams, gen] = setupGasSupply(gen, inputparams, params)
        % - inputparams is instance of ProtonicMembraneGasSupplyInputParams
        % - The setup is taken from setupInputParams method in GasSupplyGridGenerator except for the grid setup.
        % - params has the same fields as the one sent to setupInputParams method.
            
            gsgen = gen.gasSupplyGridGenerator;
            
            inputparams = gsgen.setupExternalCoupling(inputparams, params);
            
        end


        
    end
    
    
end
