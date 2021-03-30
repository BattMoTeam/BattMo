classdef LithiumBatteryInputParams < BatteryInputParams
    
    methods
        
        function paramobj = LithiumBatteryInputParams()
            paramobj = paramobj@BatteryInputParams();
            paramobj.ne.eac.am = GraphiteInputParams();
            paramobj.pe.eac.am = NMC111InputParams();
            paramobj.elyte = orgLiPF6InputParams();
        end

    end
    
end
