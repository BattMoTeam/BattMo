function paramobj = setupProtonicMembraneCellGrid(paramobj, jsonstruct)

    an    = 'Anode';
    ct    = 'Cathode';
    elyte = 'Electrolyte';

    switch jsonstruct.Geometry.type

      case '1D'

        gen  = PEMgridGenerator1D();
        
        gen.xlength  = convertUnitBattMo(jsonstruct.(elyte).length);
        gen.faceArea = jsonstruct.(elyte).faceArea;
        gen.N        = jsonstruct.(elyte).N;

        paramobj = gen.updatePEMinputParams(paramobj);
        
      otherwise
        
        error('Geometry case not recognized');
        
    end
    

    
end
