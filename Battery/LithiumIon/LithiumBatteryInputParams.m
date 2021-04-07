classdef LithiumBatteryInputParams < BatteryInputParams
    
    methods
        
        function paramobj = LithiumBatteryInputParams()
            
            paramobj = paramobj@BatteryInputParams();

            paramobj.ne = LithiumElectrodeInputParams();
            paramobj.pe = LithiumElectrodeInputParams();
            
            paramobj.ne.eac.am = GraphiteInputParams();
            paramobj.pe.eac.am = NMC111InputParams();
            
            paramobj.elyte = orgLiPF6InputParams();
            
            paramobj.SOC  = 0.5;
            paramobj.T    = 298.15;
            paramobj.J    = 0.001;
            paramobj.Ucut = 2;
            
        end

    end
    
end
