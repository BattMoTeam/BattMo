function [val, jac] = mygetfield(state, fds)
    jac = [];
    if iscell(fds)
        if numel(fds) == 1
            val = state.(fds{1});
            if iscell(val)
                val = vertcat(val{:});
            end
            if isa(val, 'ADI')
                isad = true;
                val = combineEquations(val);
                jac = val.jac{1};
                val = val.val;
            end
        else
            [val, jac] = mygetfield(state.(fds{1}), {fds{2 : end}} );
        end
    else
        error
    end
end
