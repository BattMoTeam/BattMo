classdef orgLiPF6InputParams < ElectrolyteInputParams
    
    methods

        function paramobj = orgLiPF6InputParams()
            
            paramobj = paramobj@ElectrolyteInputParams();
            
            paramobj.name = 'orgLiPF6';
            
            paramobj.compnames = {'Li', 'PF6'};
            paramobj.ncomp = 2;
            
            paramobj.indchargecarrier = 1;
            
            % Set constant values
            [~, ind] = ismember('Li', paramobj.compnames);
            tLi = 0.399;
            sp.t{ind} = tLi; % Li+ transference number, [-]
            sp.z{ind} = 1;
            
            [~, ind] = ismember('PF6', paramobj.compnames);
            sp.t{ind} = 1 - tLi; % Li+ transference number, [-]
            sp.z{ind} = -1;
            
            paramobj.sp = sp;
            
            paramobj.chargeCarrierName         = 'Li';
            paramobj.chargeCarrierFluxName     = 'LiFlux';
            paramobj.chargeCarrierSourceName   = 'LiSource';
            paramobj.chargeCarrierMassConsName = 'massCons';
            paramobj.chargeCarrierAccumName    = 'LiAccum';
            
            paramobj.sep = celgard2500InputParams();
        end
        
    end
    
end
