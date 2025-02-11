function [inputparams, gridGenerator] = setupProtonicMembraneCellGrid(inputparams, jsonstruct)

    switch jsonstruct.Geometry.type

      case '2D'

        gen = CellGridGenerator2D();
        
        gen.ly = jsonstruct.Geometry.ylength;
        gen.ny = jsonstruct.Geometry.Ny;

        gen.nxElectrolyser = jsonstruct.Electrolyser.Electrolyte.Nx;
        gen.lxElectrolyser = jsonstruct.Electrolyser.Electrolyte.xlength;

        gen.nxGasSupply = jsonstruct.GasSupply.Nx;
        gen.lxGasSupply = jsonstruct.GasSupply.xlength;
        
        inputparams = gen.updateInputParams(inputparams);
        
      otherwise

        error('Geometry case not recognized');
        
    end
    
    gridGenerator = gen;
    
end
