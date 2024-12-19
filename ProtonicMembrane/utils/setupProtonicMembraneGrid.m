function [inputparams, gridGenerator] = setupProtonicMembraneGrid(inputparams, jsonstruct)

    an    = 'Anode';
    ct    = 'Cathode';
    elyte = 'Electrolyte';

    switch jsonstruct.Geometry.type

      case '1D'

        gen  = PEMgridGenerator1D();
        
        gen.xlength  = jsonstruct.(elyte).xlength;
        gen.faceArea = jsonstruct.(elyte).faceArea;
        gen.N        = jsonstruct.(elyte).Nx;

        inputparams = gen.updateInputParams(inputparams);

      case '2D'

        gen  = PEMgridGenerator2D();
        
        gen.ylength = jsonstruct.Geometry.ylength;
        gen.Ny      = jsonstruct.Geometry.Ny;
        gen.xlength = jsonstruct.(elyte).xlength;
        gen.Nx      = jsonstruct.(elyte).Nx;
        
        inputparams = gen.updateInputParams(inputparams);
        
      otherwise

        error('Geometry case not recognized');
        
    end
    
    gridGenerator = gen;
    
end
