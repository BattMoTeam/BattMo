function [paramobj, gridGenerator] = setupBatteryGridFromJson(paramobj, jsonstruct)

    switch jsonstruct.Geometry.case

      case '1D'

        gen = BatteryGenerator1D();
        
        xlength = gen.xlength;
        
        xlength(2) = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        xlength(3) = jsonstruct.Electrolyte.Separator.thickness;
        xlength(4) = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;

        if paramobj.NegativeElectrode.include_current_collectors
            xlength(1) = jsonstruct.NegativeElectrode.CurrentCollector.thickness;
        end
        if paramobj.PositiveElectrode.include_current_collectors
            xlength(5) = jsonstruct.PositiveElectrode.CurrentCollector.thickness;
        end

        gen.xlength = xlength;
        
        gen.sepnx  = jsonstruct.Electrolyte.Separator.N;
        gen.nenx   = jsonstruct.NegativeElectrode.ActiveMaterial.N;
        gen.penx   = jsonstruct.PositiveElectrode.ActiveMaterial.N;

        if paramobj.NegativeElectrode.include_current_collectors
            gen.ccnenx = jsonstruct.NegativeElectrode.CurrentCollector.N;
        end
        if paramobj.PositiveElectrode.include_current_collectors
            gen.ccpenx = jsonstruct.PositiveElectrode.CurrentCollector.N;
        end
        
        if isfield(jsonstruct.Geometry, 'faceArea')
            gen.faceArea = jsonstruct.Geometry.faceArea;
        end

        % Now, we update the paramobj with the properties of the mesh. 
        [paramobj, gen] = gen.updateBatteryInputParams(paramobj);

      case '2D-demo'
        
      case '3D-demo'
        
        gen = BatteryGenerator3D();
        
        % zlength = gen.zlength;
        % zlength(1) = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        % zlength(2) = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        % zlength(3) = jsonstruct.Electrolyte.Separator.thickness;
        % zlength(4) = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
        % zlength(5) = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;

        % gen.sep_nz   = jsonstruct.Electrolyte.Separator.N;
        % gen.ne_am_nz = jsonstruct.NegativeElectrode.ActiveMaterial.N;
        % gen.pe_am_nz = jsonstruct.PositiveElectrode.ActiveMaterial.N;
        % gen.ne_cc_nz = jsonstruct.NegativeElectrode.CurrentCollector.N;
        % gen.pe_cc_nz = jsonstruct.PositiveElectrode.CurrentCollector.N;
        
        % Now, we update the paramobj with the properties of the mesh. 
        [paramobj, gen] = gen.updateBatteryInputParams(paramobj);
        
      case {'jellyRoll', 'sectorModel'}

        % Prepare input for SpiralBatteryGenerator.updateBatteryInputParams using json input
        
        widthDict = containers.Map();
        widthDict('ElectrolyteSeparator')     = jsonstruct.Electrolyte.Separator.thickness;
        widthDict('NegativeActiveMaterial')   = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        widthDict('NegativeCurrentCollector') = jsonstruct.NegativeElectrode.CurrentCollector.thickness;
        widthDict('PositiveActiveMaterial')   = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
        widthDict('PositiveCurrentCollector') = jsonstruct.PositiveElectrode.CurrentCollector.thickness;
        
        nwidths = [widthDict('PositiveActiveMaterial')  ; ...
                   widthDict('PositiveCurrentCollector'); ...
                   widthDict('PositiveActiveMaterial')  ; ...
                   widthDict('ElectrolyteSeparator')    ; ...
                   widthDict('NegativeActiveMaterial')  ; ...
                   widthDict('NegativeCurrentCollector'); ...
                   widthDict('NegativeActiveMaterial')  ; ...
                   widthDict('ElectrolyteSeparator')]; 
        dr = sum(nwidths);
        
        rOuter = jsonstruct.Geometry.rOuter;
        rInner = jsonstruct.Geometry.rInner;
        L      = jsonstruct.Geometry.L;
        nL     = jsonstruct.Geometry.nL;
        nas    = jsonstruct.Geometry.nas;

        nrDict = containers.Map();
        nrDict('ElectrolyteSeparator')     = jsonstruct.Electrolyte.Separator.N;
        nrDict('NegativeActiveMaterial')   = jsonstruct.NegativeElectrode.ActiveMaterial.N;
        nrDict('NegativeCurrentCollector') = jsonstruct.NegativeElectrode.CurrentCollector.N;
        nrDict('PositiveActiveMaterial')   = jsonstruct.PositiveElectrode.ActiveMaterial.N;
        nrDict('PositiveCurrentCollector') = jsonstruct.PositiveElectrode.CurrentCollector.N;

        % compute numbers of winding (this is input for spiralGrid) from outer and inner radius
        nwidths = [widthDict('PositiveActiveMaterial')  ; ...
                   widthDict('PositiveCurrentCollector'); ...
                   widthDict('PositiveActiveMaterial')  ; ...
                   widthDict('ElectrolyteSeparator')    ; ...
                   widthDict('NegativeActiveMaterial')  ; ...
                   widthDict('NegativeCurrentCollector'); ...
                   widthDict('NegativeActiveMaterial')  ; ...
                   widthDict('ElectrolyteSeparator')]; 
        dr = sum(nwidths);
        dR = rOuter - rInner; 
        % Computed number of windings
        nwindings = ceil(dR/dr);

        ne = 'NegativeElectrode';
        pe = 'PositiveElectrode';
        cc = 'CurrentCollector';
        
        params = struct('nwindings', nwindings, ...
                        'rInner'   , rInner   , ...
                        'widthDict', widthDict, ...
                        'nrDict'   , nrDict   , ...
                        'nas'      , nas      , ...
                        'L'        , L        , ...
                        'nL'       , nL       ); 

        switch jsonstruct.Geometry.case

          case 'jellyRoll'

            tabparams.NegativeElectrode = jsonstruct.(ne).(cc).tabparams;
            tabparams.PositiveElectrode = jsonstruct.(pe).(cc).tabparams;
            params.tabparams = tabparams;
            params.angleuniform = true;

            gen = SpiralBatteryGenerator();
            
          case 'sectorModel'

            gen = SectorBatteryGenerator();

          otherwise

            error('Geometry.case not recognized');
            
        end
        
        [paramobj, gen] = gen.updateBatteryInputParams(paramobj, params);
        
      otherwise
        
        error('Geometry case not recognized')
        
    end

    gridGenerator = gen;

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
