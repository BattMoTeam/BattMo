function jsonstruct = setJsonStructField(jsonstruct, fieldnamelist, value, varargin)

    if ischar(fieldnamelist)
        % handle case wher fieldnamelist is just a char
        fieldnamelist = {fieldnamelist};
        jsonstruct = setJsonStructField(jsonstruct, fieldnamelist, value, varargin{:});
        return
    end

    fieldname = fieldnamelist{1};

    if numel(fieldnamelist) > 1
        fieldnamelist = fieldnamelist(2 : end);
        setValue = false;
    else
        setValue = true;
    end

    if ~isempty(jsonstruct) && isfield(jsonstruct, fieldname)

        if setValue

            currentValue = jsonstruct.(fieldname);

            equalValue = isequal(currentValue, value);
            
            if equalValue
                
                % do nothing. Value was already set to given value
                return
                
            else
                
                opt = struct('handleMisMatch', 'error');
                opt = merge_options(opt, varargin{:});
                
                switch opt.handleMisMatch

                  case 'quiet'

                    jsonstruct.(fieldname) = value;
                    
                  case 'warning'

                    fprintf('mismatch values in assignment of %s. We use the given value\n', fieldname)
                    jsonstruct.(fieldname) = value;
                    
                  case 'error'

                    errortxt = sprintf('mismatch values in assignment of %s. We use the given value\n', fieldname);
                    error(errortxt);
                    
                end

            end
        else

            jsonstruct.(fieldname) = setJsonStructField(jsonstruct.(fieldname), fieldnamelist, value, varargin{:});

        end

    else

        if setValue

            jsonstruct.(fieldname) = value;

        else

            jsonstruct.(fieldname) = setJsonStructField([], fieldnamelist, value, varargin{:});
            
        end
        
    end
        
    
end