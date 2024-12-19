classdef CellGridGenerator

    properties

        parentGrid
        
        % Use to generate the structure of the internal coupling terms
        electrolyserGridGenerator
        gasSupplyGridGenerator
        
    end

    methods
        
        function [inputparams, gen] = updateInputParams(gen, inputparams, params)
            
            error('virtual function');
            
        end
        
        function [inputparams, gen] = setupInputParams(gen, inputparams, params)
        % This is the main function which should be called by the derived class at the end of the :meth:`updateInputParams` method

            
            [inputparams, gen] = gen.setupGrid(inputparams, params);

            gen = gen.setupElectrolyserGridGenerator(inputparams, params);
            gen = gen.setupGasSupplyGridGenerator(inputparams, params);

            inputparams.Electrolyser.G = genSubGrid(gen.parentGrid, params.Electrolyser.cellinds);
            inputparams.GasSupply.G    = genSubGrid(gen.parentGrid, params.GasSupply.cellinds);
            
            [inputparams.Electrolyser, gen] = gen.setupElectrolyser(inputparams.Electrolyser, params.Electrolyser);
            [inputparams.GasSupply, gen]    = gen.setupGasSupply(inputparams.GasSupply, params.GasSupply);

            [inputparams, gen] = gen.setupElectrolyserGasSupplyCoupling(inputparams);
            
        end

        function [inputparams, gen] = setupGrid(gen, inputparams, params)
            
            error('virtual function');
            
        end
        
        function gen = setupElectrolyserGridGenerator(gen, inputparams, params)
        % setup gen.electrolyserGridGenerator

            error('virtual function');
            
        end

        function gen = setupGasSupplyGridGenerator(gen, inputparams, params)
        % setup gen.gasSupplyGridGenerator
            
            error('virtual function');
            
        end

        function [inputparams, gen] = setupElectrolyser(gen, inputparams, params)
        % - inputparams is instance of ProtonicMembraneInputParams
        % - The setup is taken from setupInputParams method in PEMgridGenerator except for the grid setup.
        %   Note that we assume that electrolyte and PEM have same grid (this may changed...)
        % - params has the same fields as the one sent to setupInputParams method.
        %   params.Electrolyte.cellind should correspond to gen.parentGrid (which is same as gen.electrolyserGridGenerator.parentGrid).

            cgen = gen.electrolyserGridGenerator;
            
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
