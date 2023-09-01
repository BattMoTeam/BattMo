function state=getStateFromVariables(y,model,dims,state0)
    varNames=model.primaryVariableNames;
    state=state0;
    nbegin=1;
    for i = 1:length(varNames)
        nend=nbegin + dims(i)-1;

%Add fields to EMPTY struct
%         add=addFields(varNames{i},y(nbegin:nend));
% 
%         if ~isfield(state,varNames{i}{1})
%             state.(varNames{i}{1})=add;
%         else
%             k=2;
%             tmp=state.(varNames{i}{1});
%             while k<length(varNames{i}) && isfield(tmp,varNames{i}{k})
%                 tmp=tmp.(varNames{i}{k});
%                 k=k+1;
%             end
%             tmp.(varNames{i}{k})=y(nbegin:nend);
% 
%             for j=k-1:-1:2
%                 tmp=struct(varNames{i}{j},tmp);
%             end
%             state.(varNames{i}{1})=tmp;
%         end

        call='state';
        for j=1:length(varNames{i})
            call=strcat(call,'.',varNames{i}{j});
        end
        call=strcat(call,"=y(nbegin:nend);");
        eval(call);
        nbegin=nend+1;
    end
%     state.time=state0.time;
%     state.ThermalModel=state0.ThermalModel;
    state=model.addVariables(state);
end

