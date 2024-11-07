function jsonstruct = equalizeJsonStructFields(jsonstruct, fieldnamelists, varargin)
% same as equalizeJsonStructField (without 's' at end) but iterate over the fieldnamelists;

    fielnamelist1 = fieldnamelists{1};

    for ifl = 2 : numel(fieldnamelists)
        
        fieldnamlist2 = fieldnamelists{ifl};

        jsonstruct = equalizeJsonStructField(jsonstruct, fielnamelist1, fieldnamlist2, varargin{:});
        
    end

    
end
