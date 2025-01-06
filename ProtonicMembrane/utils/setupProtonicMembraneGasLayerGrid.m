function [inputparams, gridGenerator] = setupProtonicMembraneGasLayerGrid(inputparams, jsonstruct)

    geom = 'Geometry';
    
    switch jsonstruct.Geometry.type

      case '2D'

        gen = GasSupplyGridGenerator2D();

        gen.Nx      = jsonstruct.(geom).Nx;
        gen.Ny      = jsonstruct.(geom).Ny;
        gen.xlength = jsonstruct.(geom).xlength;
        gen.ylength = jsonstruct.(geom).ylength;
        
      case '1D'
        
        gen = GasSupplyGridGenerator1D();

        gen.Nx       = jsonstruct.(geom).Nx;
        gen.xlength  = jsonstruct.(geom).xlength;
        gen.faceArea = jsonstruct.(geom).faceArea;
        
      otherwise

        error('dimcase not recognized')
        
    end

    inputparams = gen.updateInputParams(inputparams);

    gridGenerator = gen;
    
end
