function [inputparams, gridGenerator] = setupBatteryGridFromJson(inputparams, jsonstruct)

    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    elyte = 'Electrolyte';
    co    = 'Coating';
    sep   = 'Separator';
    cc    = 'CurrentCollector';

    switch jsonstruct.Geometry.case

      case '1D'

        gen = BatteryGeneratorP2D();

        xlength = gen.xlength;

        xlength(2) = jsonstruct.NegativeElectrode.Coating.thickness;
        xlength(3) = jsonstruct.Separator.thickness;
        xlength(4) = jsonstruct.PositiveElectrode.Coating.thickness;

        if inputparams.NegativeElectrode.include_current_collectors
            xlength(1) = jsonstruct.NegativeElectrode.CurrentCollector.thickness;
        end
        if inputparams.PositiveElectrode.include_current_collectors
            xlength(5) = jsonstruct.PositiveElectrode.CurrentCollector.thickness;
        end

        gen.xlength = xlength;

        gen.sepnx  = jsonstruct.Separator.N;
        gen.nenx   = jsonstruct.NegativeElectrode.Coating.N;
        gen.penx   = jsonstruct.PositiveElectrode.Coating.N;

        if inputparams.NegativeElectrode.include_current_collectors
            gen.ccnenx = jsonstruct.NegativeElectrode.CurrentCollector.N;
        end
        if inputparams.PositiveElectrode.include_current_collectors
            gen.ccpenx = jsonstruct.PositiveElectrode.CurrentCollector.N;
        end

        if isfield(jsonstruct.Geometry, 'faceArea')
            gen.faceArea = jsonstruct.Geometry.faceArea;
        end

        if isfield(jsonstruct.Geometry, 'resolutionFactor')
            gen.resolutionFactor = jsonstruct.Geometry.resolutionFactor;
            gen = gen.applyResolutionFactors();
        end

        % Now, we update the inputparams with the properties of the grid.
        [inputparams, gen] = gen.updateBatteryInputParams(inputparams);

      case 'multiLayerPouch'

        geom = 'Geometry';

        gen = BatteryGeneratorMultilayerPouch();

        gen.pouch_width  = jsonstruct.(geom).width;
        gen.pouch_height = jsonstruct.(geom).length;

        gen.unit_cell_thickness = [jsonstruct.(ne).(cc).thickness ; ...
                                   jsonstruct.(ne).(co).thickness ; ...
                                   jsonstruct.(sep).thickness     ; ...
                                   jsonstruct.(pe).(co).thickness ; ...
                                   jsonstruct.(pe).(cc).thickness];

        gen.sep_nz   = jsonstruct.(sep).N;
        gen.ne_co_nz = jsonstruct.(ne).(co).N;
        gen.pe_co_nz = jsonstruct.(pe).(co).N;
        gen.ne_cc_nz = jsonstruct.(ne).(cc).N;
        gen.pe_cc_nz = jsonstruct.(pe).(cc).N;


        gen.tab_width     = jsonstruct.(geom).tab.width;
        gen.ne_tab_height = jsonstruct.(geom).tab.(ne).length;
        gen.pe_tab_height = jsonstruct.(geom).tab.(pe).length;

        gen.n_layers = jsonstruct.(geom).nLayers;

        gen.elyte_nx  = jsonstruct.(geom).(elyte).Nx;
        gen.elyte_ny  = jsonstruct.(geom).(elyte).Ny;
        gen.tab_nx    = jsonstruct.(geom).tab.Nx;
        gen.ne_tab_ny = jsonstruct.(geom).tab.(ne).Ny;
        gen.pe_tab_ny = jsonstruct.(geom).tab.(pe).Ny;

        if isfield(jsonstruct.(geom).tab, 'cap_tabs') && jsonstruct.(geom).tab.cap_tabs
            gen.cap_tabs = true;
        else
            gen.cap_tabs = false;
        end
        % Now, we update the inputparams with the properties of the grid.
        [inputparams, gen] = gen.updateBatteryInputParams(inputparams);

      case '2D-demo'

      case '3D-demo'

        gen = BatteryGeneratorP4D();

        zlength = gen.zlength;
        zlength(1) = jsonstruct.NegativeElectrode.Coating.thickness;
        zlength(2) = jsonstruct.NegativeElectrode.Coating.thickness;
        zlength(3) = jsonstruct.Separator.thickness;
        zlength(4) = jsonstruct.PositiveElectrode.Coating.thickness;
        zlength(5) = jsonstruct.PositiveElectrode.Coating.thickness;
        gen.zlength = zlength;

        gen.sep_nz   = jsonstruct.Separator.N;
        gen.ne_co_nz = jsonstruct.NegativeElectrode.Coating.N;
        gen.pe_co_nz = jsonstruct.PositiveElectrode.Coating.N;
        gen.ne_cc_nz = jsonstruct.NegativeElectrode.CurrentCollector.N;
        gen.pe_cc_nz = jsonstruct.PositiveElectrode.CurrentCollector.N;

        xlength = gen.xlength;
        xlength(1) = jsonstruct.NegativeElectrode.CurrentCollector.tab.width;
        xlength(3) = jsonstruct.PositiveElectrode.CurrentCollector.tab.width;
        xlength(2) = jsonstruct.Geometry.width - (xlength(1) + xlength(3));
        gen.xlength = xlength;

        gen.ne_cc_nx     = jsonstruct.NegativeElectrode.CurrentCollector.tab.Nw;
        gen.pe_cc_nx     = jsonstruct.PositiveElectrode.CurrentCollector.tab.Nw;
        gen.int_elyte_nx = jsonstruct.Geometry.Nw - (gen.ne_cc_nx + gen.pe_cc_nx);

        ylength = gen.ylength;
        ylength(1) = jsonstruct.NegativeElectrode.CurrentCollector.tab.height;
        ylength(3) = jsonstruct.PositiveElectrode.CurrentCollector.tab.height;
        ylength(2) = jsonstruct.Geometry.height;
        gen.ylength = ylength;

        gen.ne_cc_ny = jsonstruct.NegativeElectrode.CurrentCollector.tab.Nh;
        gen.pe_cc_ny = jsonstruct.PositiveElectrode.CurrentCollector.tab.Nh;
        gen.elyte_ny = jsonstruct.Geometry.Nh;

        gen.externalHeatTransferCoefficient = jsonstruct.ThermalModel.externalHeatTransferCoefficient;
        if isfield(jsonstruct.ThermalModel, 'externalHeatTransferCoefficientTab')
            gen.externalHeatTransferCoefficientTab = jsonstruct.ThermalModel.externalHeatTransferCoefficientTab;
        else
            gen.externalHeatTransferCoefficientTab = gen.externalHeatTransferCoefficient;
        end

        % Now, we update the inputparams with the properties of the grid.
        [inputparams, gen] = gen.updateBatteryInputParams(inputparams);

      case {'jellyRoll', 'sectorModel'}

        % Prepare input for SpiralBatteryGenerator.updateBatteryInputParams using json input

        widthDict = containers.Map();
        widthDict('Separator')                = jsonstruct.Separator.thickness;
        widthDict('NegativeCoating')          = jsonstruct.NegativeElectrode.Coating.thickness;
        widthDict('NegativeCurrentCollector') = jsonstruct.NegativeElectrode.CurrentCollector.thickness;
        widthDict('PositiveCoating')          = jsonstruct.PositiveElectrode.Coating.thickness;
        widthDict('PositiveCurrentCollector') = jsonstruct.PositiveElectrode.CurrentCollector.thickness;

        nwidths = [widthDict('PositiveCoating')          ; ...
                   widthDict('PositiveCurrentCollector') ; ...
                   widthDict('PositiveCoating')          ; ...
                   widthDict('Separator')                ; ...
                   widthDict('NegativeCoating')          ; ...
                   widthDict('NegativeCurrentCollector') ; ...
                   widthDict('NegativeCoating')          ; ...
                   widthDict('Separator')];
        dr = sum(nwidths);

        rOuter = jsonstruct.Geometry.rOuter;
        rInner = jsonstruct.Geometry.rInner;
        L      = jsonstruct.Geometry.L;
        nL     = jsonstruct.Geometry.nL;
        nas    = jsonstruct.Geometry.nas;

        if isfield(jsonstruct.Geometry, 'refLcoef')
            refLcoef = jsonstruct.Geometry.refLcoef;
        else
            refLcoef = [];
        end

        if isfield(jsonstruct.Geometry, 'exteriorNegativeElectrodeLayer')
            exteriorNegativeElectrodeLayer = jsonstruct.Geometry.exteriorNegativeElectrodeLayer;
        else
            exteriorNegativeElectrodeLayer = false;
        end
        
        nrDict = containers.Map();
        nrDict('Separator')                = jsonstruct.Separator.N;
        nrDict('NegativeCoating')          = jsonstruct.NegativeElectrode.Coating.N;
        nrDict('NegativeCurrentCollector') = jsonstruct.NegativeElectrode.CurrentCollector.N;
        nrDict('PositiveCoating')          = jsonstruct.PositiveElectrode.Coating.N;
        nrDict('PositiveCurrentCollector') = jsonstruct.PositiveElectrode.CurrentCollector.N;

        % compute numbers of winding (this is input for spiralGrid) from outer and inner radius
        nwidths = [widthDict('PositiveCoating')          ; ...
                   widthDict('PositiveCurrentCollector') ; ...
                   widthDict('PositiveCoating')          ; ...
                   widthDict('Separator')                ; ...
                   widthDict('NegativeCoating')          ; ...
                   widthDict('NegativeCurrentCollector') ; ...
                   widthDict('NegativeCoating')          ; ...
                   widthDict('Separator')];
        dr = sum(nwidths);
        dR = rOuter - rInner;
        % Computed number of windings
        nwindings = ceil(dR/dr);

        params = struct('nwindings' , nwindings, ...
                        'rInner'    , rInner   , ...
                        'widthDict' , widthDict, ...
                        'nrDict'    , nrDict   , ...
                        'nas'       , nas      , ...
                        'L'         , L        , ...
                        'nL'        , nL       , ...
                        'refLcoef'  , refLcoef , ...
                        'exteriorNegativeElectrodeLayer', exteriorNegativeElectrodeLayer );

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

        [inputparams, gen] = gen.updateBatteryInputParams(inputparams, params);

      otherwise

        error('Geometry case not recognized')

    end

    gridGenerator = gen;

end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
