function jsonstruct = setJsonStructField(jsonstruct, fieldnamelist, value, varargin)
    % in varargin, the key 'handleMisMatch' can take following vaules
    % - 'error'   : returns error if value already set and does not match with the new given one (default)
    % - 'quiet'   : does not warn about case above
    % - 'warning' : warn but no error

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

    if ~isempty(jsonstruct) && (isfield(jsonstruct, fieldname) || isprop(jsonstruct, fieldname))

        if setValue

            currentValue = jsonstruct.(fieldname);

            equalValue = isequal(currentValue, value);
            
            if equalValue
                
                % do nothing. Value was already set to given value
                return
                
            else
                
                opt = struct('handleMisMatch', 'error', ...
                             'errorMessage', []);
                opt = merge_options(opt, varargin{:});
                
                switch opt.handleMisMatch

                  case 'quiet'

                    jsonstruct.(fieldname) = value;
                    
                  case 'warning'

                    fprintf('mismatch values in assignment of %s. We use the given value\n', fieldname)
                    jsonstruct.(fieldname) = value;
                    
                  case 'error'

                    if isempty(opt.errorMessage)
                        errorMessage = sprintf('mismatch values in assignment of %s. We use the given value\n', fieldname);
                    else
                        errorMessage = opt.errorMessage;
                    end
                    error(errorMessage);
                    
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
