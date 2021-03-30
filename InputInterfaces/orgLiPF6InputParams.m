classdef orgLiPF6InputParams < ElectrolyteInputParams
    
    methods

        function paramobj = orgLiPF6InputParams()
            
            paramobj = paramobj@ElectrolyteInputParams();
            
            paramobj.compnames = {'Li', 'PF6'};
            paramobj.ncomp = 2;
            
            % Set constant values
            [~, ind] = ismember('Li', paramobj.compnames);
            tLi = 0.399;
            sp.t{ind} = tLi; % Li+ transference number, [-]
            sp.z{ind} = 1;
            
            [~, ind] = ismember('PF6', paramobj.compnames);
            sp.t{ind} = 1 - tLi; % Li+ transference number, [-]
            sp.z{ind} = -1;
            
            paramobj.sp = sp;
            
        end
        
    end
    
end
