classdef GasSupplyPEMgridGenerator

    properties

        G
        
        % Use to generate the structure of the internal coupling terms
        cellGridGenerator
        gasSupplyGridGenerator
        
    end

    methods
        
        function [paramobj, gen] = updateInputParams(gen, paramobj, params)
            
            error('virtual function');
            
        end
        
        function [paramobj, gen] = setupInputParams(gen, paramobj, params)
        % This is the main function which should be called by the derived class at the end of the :meth:`updateInputParams` method

            
            [paramobj, gen] = gen.setupGrid(paramobj, params);

            [paramobj, gen] = gen.setupCellGridGenerator(paramobj, params);
            [paramobj, gen] = gen.setupGasSupplyGridGenerator(paramobj, params);

            [paramobj.Cell, gen]      = gen.setupCell(paramobj.Cell, params.Cell);
            [paramobj.GasSupply, gen] = gen.setupGasSupply(paramobj.GasSupply, params.GasSupply);

            [paramobj, gen] = gen.setupCellGasSupplyCoupling(paramobj);
            
        end

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
            
            error('virtual function');
            
        end
        
        function [paramobj, gen] = setupCellGridGenerator(gen, paramobj, params)
        % setup cellGridGenerator

            error('virtual function');
            
        end

        function [paramobj, gen] = setupGasSupplyGridGenerator(gen, paramobj, params)
        % setup gasSupplyGridGenerator
            
            error('virtual function');
            
        end

        function [paramobj, gen] = setupCellGrid(gen, paramobj, params)
        % paramobj belongs to CellGeneratorInputParams

            globG    = gen.G;
            cellinds = params.cellinds;

            G = genSubGrid(globG, cellinds);

            paramobj.G = G;
            gen.cellGridGenerator.G = G;
            
        end

        function [paramobj, gen] = setupGasSupplyGrid(gen, paramobj, params)
        % paramobj belongs to GasSupplyGridGeneratorInputParams
            
            globG    = gen.G;
            cellinds = params.cellinds;

            G = genSubGrid(globG, cellinds);

            paramobj.G = G;
            gen.gasSupplyGridGenerator.G = G;
            
        end

        function [paramobj, gen] = setupCell(gen, paramobj, params)

            [paramobj, gen] = setupCellGrid(gen, paramobj, params);
            
            cgen = gen.cellGridGenerator;
            params_elyte = pickField(params, 'Electrolyte');
            paramobj.Electrolyte = cgen.setupElectrolyte(paramobj.Electrolyte, params_elyte);
            paramobj = cgen.setupElectrodeElectrolyteCoupTerm(paramobj, params);

        end

        function [paramobj, gen] = setupGasSupply(gen, paramobj, params)
        % paramobj belongs to GasSupplyInputParams
            
            [paramobj, gen] = setupGasSupplyGrid(gen, paramobj, params);
            
            glgen = gen.gasSupplyGridGenerator;
            [paramobj, glgen] = glgen.setupExternalCoupling(paramobj, params);
            
            gen.gasSupplyGridGenerator = glgen;
            
        end
        
    end
    
    
end
