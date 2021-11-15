classdef Separator < BaseModel
    
    properties
        
        % Physicochemical properties
        thickness      % Thickness,        [m]
        porosity       % Porosity,         [-]
        rp             % Pore radius,      [m]
        Gurley         % Gurley number,    [s]
        volumeFraction % Volume fraction,  [-]
        
        thermalConductivity % intrinsic thermal conductivity value
        heatCapacity % intrinsic heat capacity value

        EffectiveThermalConductivity
        EffectiveHeatCapacity
        density % [kg m^-3]        
    end
    
    methods
        function model = Separator(paramobj)
            
            model = model@BaseModel();

            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            % in the case of the separator, probably this does not matter as no computation is actually done on this grid
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G'                  , ...
                       'thickness'          , ...
                       'porosity'           , ...
                       'rp'                 , ...
                       'Gurley'             , ...
                       'thermalConductivity', ...
                       'heatCapacity'       , ...
                       'density'};
            model = dispatchParams(model, paramobj, fdnames);
            
            model.volumeFraction = 1 - model.porosity;
            model.EffectiveThermalConductivity = model.thermalConductivity.*(model.volumeFraction).^1.5;
            model.EffectiveHeatCapacity = model.heatCapacity.*model.volumeFraction;
            
        end
    end
    
end

