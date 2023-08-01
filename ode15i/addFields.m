% function ret=addFields(names,obj,final)
%     if length(names)>1 && ~ischar(names)
%         tmp=struct(names{1},obj);
%         ret=addFields({names{2:end}},tmp.(names{1}),final);
%     else
%         obj.(names{1})=final;
%         ret=obj;
%     end
% end

function ret=addFields(names,obj)
    if length(names)>1 && ~ischar(names)
        ret=addFields({names{1:end-1}},struct(names{end},obj));
    else
        ret=obj;
    end
end