function paramobj = setupElectrolyserGridFromJson(paramobj, jsonstruct)

    inm = 'IonomerMembrane';
    her = 'HydrogenEvolutionElectrode';
    oer = 'OxygenEvolutionElectrode';
    ptl = 'PorousTransportLayer';
    ctl = 'CatalystLayer';
    
    switch jsonstruct.Geometry.case

      case '1D'

        xs = nan(5, 1);
        
        xs(1) = jsonstruct.(oer).(ptl).length;
        xs(2) = jsonstruct.(oer).(ctl).length;
        xs(3) = jsonstruct.(inm).length;
        xs(4) = jsonstruct.(her).(ctl).length;        
        xs(5) = jsonstruct.(her).(ptl).length;

        cxs = [jsonstruct.(oer).(ptl).cellsize;
               jsonstruct.(oer).(ctl).cellsize;
               jsonstruct.(inm).cellsize;
               jsonstruct.(her).(ctl).cellsize;        
               jsonstruct.(her).(ptl).cellsize];

       gen = ElectrolyserGridGenerator1D(xs, cxs);
       paramobj = gen.updateElectrolyserInputParams(paramobj);

       % We add the solid volume fraction 

       vf_inm = ones(paramobj.(inm).G.cells.num, 1);

       eldes = {oer, her};
       
       for ielde = 1 : numel(eldes)
           
           elde = eldes{ielde};
           svf_elde = jsonstruct.(elde).(ptl).widthFraction;
           svf_elde = svf_elde*ones(paramobj.(elde).(ptl).G.cells.num, 1);
           
           coupterm = paramobj.(elde).couplingTerm;
           w_elde = jsonstruct.(elde).(ctl).widthFraction;
           switch elde
             case oer
               side = 'oxygenCatalystSide';
             case her
               side = 'hydrogenCatalystSide';
           end
           w_inm = jsonstruct.(inm).widthFractions.(side);
           svf_elde(coupterm.couplingcells(:, 1)) = w_elde + w_inm;
           
           coupterms = paramobj.couplingTerms;
           coupnames = cellfun(@(x) x.name, coupterms, 'uniformoutput', false);
           coupterm = getCoupTerm(coupterms, sprintf('%s-%s', elde, inm), coupnames);
           
           vf_inm(coupterm.couplingcells(:, 2)) = w_inm;

           paramobj.(elde).(ptl).solidVolumeFraction = svf_elde;
           
       end

       paramobj.(inm).volumeFraction = vf_inm;
       
      otherwise
        
        error('Geometry case not implemented');
        
    end
        
end

