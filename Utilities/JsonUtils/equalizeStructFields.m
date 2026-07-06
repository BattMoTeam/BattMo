function jsonstruct = equalizeStructFields(jsonstruct, fieldnamelists, varargin)
% same as equalizeStructField (without 's' at end) but iterate over the fieldnamelists;

    fielnamelist1 = fieldnamelists{1};

    for ifl = 2 : numel(fieldnamelists)
        
        fieldnamlist2 = fieldnamelists{ifl};

        jsonstruct = equalizeStructField(jsonstruct, fielnamelist1, fieldnamlist2, varargin{:});
        
    end

    
end
