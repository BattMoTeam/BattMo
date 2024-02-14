function [inputparams, gridGenerator] = setupProtonicMembraneCellGrid(inputparams, jsonstruct)

    an    = 'Anode';
    ct    = 'Cathode';
    elyte = 'Electrolyte';

    switch jsonstruct.Geometry.type

      case '1D'

        gen  = PEMgridGenerator1D();
        
        gen.xlength  = convertUnitBattMo(jsonstruct.(elyte).xlength);
        gen.faceArea = jsonstruct.(elyte).faceArea;
        gen.N        = jsonstruct.(elyte).Nx;

        inputparams = gen.updatePEMinputParams(inputparams);

      case '2D'

        gen  = PEMgridGenerator2D();
        
        gen.ylength = convertUnitBattMo(jsonstruct.Geometry.ylength);
        gen.Ny      = jsonstruct.Geometry.Ny;
        gen.xlength = convertUnitBattMo(jsonstruct.(elyte).xlength);
        gen.Nx      = jsonstruct.(elyte).Nx;
        
        inputparams = gen.updatePEMinputParams(inputparams);
        
      otherwise

        error('Geometry case not recognized');
        
    end
    
    gridGenerator = gen;
    
end
