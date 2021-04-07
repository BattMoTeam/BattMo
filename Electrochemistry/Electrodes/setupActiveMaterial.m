function am = setupActiveMaterial(paramobj)

    amName = paramobj.name;
    
    switch amName 
      case 'NMC111'
        am = NMC111(paramobj);
      case 'Graphite'
        am = Graphite(paramobj);
      otherwise
        error('active material not recognized');
    end
    
end
