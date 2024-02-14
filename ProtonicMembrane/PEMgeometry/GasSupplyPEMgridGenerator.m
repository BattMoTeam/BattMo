classdef GasSupplyPEMgridGenerator

    properties

        G
        
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

            [inputparams, gen] = gen.setupCellGridGenerator(inputparams, params);
            [inputparams, gen] = gen.setupGasSupplyGridGenerator(inputparams, params);

            [inputparams.Cell, gen]      = gen.setupCell(inputparams.Cell, params.Cell);
            [inputparams.GasSupply, gen] = gen.setupGasSupply(inputparams.GasSupply, params.GasSupply);

            [inputparams, gen] = gen.setupCellGasSupplyCoupling(inputparams);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
            
            error('virtual function');
            
        end
        
        function [inputparams, gen] = setupCellGridGenerator(gen, inputparams, params)
        % setup cellGridGenerator

            error('virtual function');
            
        end

        function [inputparams, gen] = setupGasSupplyGridGenerator(gen, inputparams, params)
        % setup gasSupplyGridGenerator
            
            error('virtual function');
            
        end

        function [inputparams, gen] = setupCellGrid(gen, inputparams, params)
        % inputparams belongs to CellGeneratorInputParams

            globG    = gen.G;
            cellinds = params.cellinds;

            G = genSubGrid(globG, cellinds);

            inputparams.G = G;
            gen.cellGridGenerator.G = G;
            
        end

        function [inputparams, gen] = setupGasSupplyGrid(gen, inputparams, params)
        % inputparams belongs to GasSupplyGridGeneratorInputParams
            
            globG    = gen.G;
            cellinds = params.cellinds;

            G = genSubGrid(globG, cellinds);

            inputparams.G = G;
            gen.gasSupplyGridGenerator.G = G;
            
        end

        function [inputparams, gen] = setupCell(gen, inputparams, params)

            [inputparams, gen] = setupCellGrid(gen, inputparams, params);
            
            cgen = gen.cellGridGenerator;
            params_elyte = pickField(params, 'Electrolyte');
            inputparams.Electrolyte = cgen.setupElectrolyte(inputparams.Electrolyte, params_elyte);
            inputparams = cgen.setupElectrodeElectrolyteCoupTerm(inputparams, params);

        end

        function [inputparams, gen] = setupGasSupply(gen, inputparams, params)
        % inputparams belongs to GasSupplyInputParams
            
            [inputparams, gen] = setupGasSupplyGrid(gen, inputparams, params);
            
            glgen = gen.gasSupplyGridGenerator;
            [inputparams, glgen] = glgen.setupExternalCoupling(inputparams, params);
            
            gen.gasSupplyGridGenerator = glgen;
            
        end
        
    end
    
    
end
