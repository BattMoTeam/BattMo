function paramobj = setupBatteryGridFromJson(paramobj, jsonstruct)

    switch jsonstruct.Geometry.case

      case '1D'

        gen = BatteryGenerator1D();
        
        xlength = gen.xlength;
        xlength(2) = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        xlength(3) = jsonstruct.Electrolyte.Separator.thickness;
        xlength(4) = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
        if paramobj.NegativeElectrode.include_current_collectors
            xlength(1) = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        end
        if paramobj.PositiveElectrode.include_current_collectors
            xlength(5) = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
        end
        if isfield(jsonstruct.Geometry, 'faceArea')
            gen.faceArea = jsonstruct.Geometry.faceArea;
        end

        gen.sepnx  = jsonstruct.NegativeElectrode.ActiveMaterial.N;
        gen.nenx   = jsonstruct.Electrolyte.Separator.N;
        gen.penx   = jsonstruct.PositiveElectrode.ActiveMaterial.N;
        
        % Now, we update the paramobj with the properties of the mesh. 
        paramobj = gen.updateBatteryInputParams(paramobj);

      case {'2D-demo', '3d-demo'}

        error('not yet implemented');

      case {'jellyRoll', 'sectorModel'}

        % Prepare input for SpiralBatteryGenerator.updateBatteryInputParams using json input
        
        widthDict = containers.Map();
        widthDict('ElectrolyteSeparator')     = jsonstruct.Electrolyte.Separator.thickness;
        widthDict('NegativeActiveMaterial')   = jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
        widthDict('NegativeCurrentCollector') = jsonstruct.NegativeElectrode.CurrentCollector.thickness;
        widthDict('PositiveActiveMaterial')   = jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
        widthDict('PositiveCurrentCollector') = jsonstruct.PositiveElectrode.CurrentCollector.thickness;
        
        nwidths = [widthDict('PositiveActiveMaterial');...
                   widthDict('PositiveCurrentCollector');...
                   widthDict('PositiveActiveMaterial');...
                   widthDict('ElectrolyteSeparator');...
                   widthDict('NegativeActiveMaterial');...
                   widthDict('NegativeCurrentCollector');...
                   widthDict('NegativeActiveMaterial');...
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
        nwidths = [widthDict('PositiveActiveMaterial');...
                   widthDict('PositiveCurrentCollector');...
                   widthDict('PositiveActiveMaterial');...
                   widthDict('ElectrolyteSeparator');...
                   widthDict('NegativeActiveMaterial');...
                   widthDict('NegativeCurrentCollector');...
                   widthDict('NegativeActiveMaterial');...
                   widthDict('ElectrolyteSeparator')]; 
        dr = sum(nwidths);dR = rOuter - rInner; 
        % Computed number of windings
        nwindings = ceil(dR/dr);

        ne = 'NegativeElectrode';
        pe = 'PositiveElectrode';
        cc = 'CurrentCollector';
        
        tabparams.NegativeElectrode = jsonstruct.(ne).(cc).tabparams;
        tabparams.PositiveElectrode = jsonstruct.(pe).(cc).tabparams;

        spiralparams = struct('nwindings'   , nwindings, ...
                              'rInner'      , rInner   , ...
                              'widthDict'   , widthDict, ...
                              'nrDict'      , nrDict   , ...
                              'nas'         , nas      , ...
                              'L'           , L        , ...
                              'nL'          , nL       , ...
                              'tabparams'   , tabparams, ...
                              'angleuniform', true); 

        switch jsonstruct.Geometry.case

          case 'jellyRoll'

            gen = SpiralBatteryGenerator();
            
          case 'sectorModel'

            gen = SectorBatteryGenerator();
            
        end
        
        paramobj = gen.updateBatteryInputParams(paramobj, spiralparams);

        
      otherwise
        
        error('Geometry case not recognized')
        
    end

end
