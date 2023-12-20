function inputparams = setupElectrolyserGridFromJson(inputparams, jsonstruct)

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
       inputparams = gen.updateElectrolyserInputParams(inputparams);

       % We add the solid volume fraction 

       vf_inm = ones(inputparams.(inm).G.cells.num, 1);

       eldes = {oer, her};
       
       for ielde = 1 : numel(eldes)
           
           elde = eldes{ielde};
           svf_elde = jsonstruct.(elde).(ptl).widthFraction;
           svf_elde = svf_elde*ones(inputparams.(elde).(ptl).G.cells.num, 1);
           
           coupterm = inputparams.(elde).couplingTerm;
           w_elde = jsonstruct.(elde).(ctl).widthFraction;
           switch elde
             case oer
               side = 'oxygenCatalystSide';
             case her
               side = 'hydrogenCatalystSide';
           end
           w_inm = jsonstruct.(inm).widthFractions.(side);
           svf_elde(coupterm.couplingcells(:, 1)) = w_elde + w_inm;
           
           coupterms = inputparams.couplingTerms;
           coupnames = cellfun(@(x) x.name, coupterms, 'uniformoutput', false);
           coupterm = getCoupTerm(coupterms, sprintf('%s-%s', elde, inm), coupnames);
           
           vf_inm(coupterm.couplingcells(:, 2)) = w_inm;

           inputparams.(elde).(ptl).solidVolumeFraction = svf_elde;
           
       end

       inputparams.(inm).volumeFraction = vf_inm;
       
      otherwise
        
        error('Geometry case not implemented');
        
    end
        
end




%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
