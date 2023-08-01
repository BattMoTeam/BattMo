function [y,dims]=getODEVariables(state,model)
    varNames=model.primaryVariableNames;
    y=[];
    dims=[];
    for i = 1:length(varNames)
        tmp=state;
        for j = 1:length(varNames{i})
            tmp=tmp.(varNames{i}{j});
        end
        dims=vertcat(dims,length(tmp.val));
        y=vertcat(y,tmp);
    end
end

