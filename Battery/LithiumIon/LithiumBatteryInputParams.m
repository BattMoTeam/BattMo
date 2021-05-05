classdef LithiumBatteryInputParams < BatteryInputParams
    
    methods
        
        function paramobj = LithiumBatteryInputParams()
            
            paramobj = paramobj@BatteryInputParams();

            paramobj.ne = LithiumElectrodeInputParams();
            paramobj.pe = LithiumElectrodeInputParams();
            
            % Once the electrode are instantiated, we fill in the active material inputs
            paramobj.ne.eac.am = GraphiteInputParams();
            paramobj.pe.eac.am = NMC111InputParams();
            
            paramobj.elyte = orgLiPF6InputParams();
            
            % some (dummy) values for energy equation
            tC = 1e-8;
            hC = 1e5;
            
            paramobj.elyte.thermalConductivity = tC;
            paramobj.elyte.heatCapacity = hC;
            paramobj.elyte.sep.thermalConductivity = tC;
            paramobj.elyte.sep.heatCapacity = hC;            
            
            eldes = {'pe', 'ne'};
            cpts = {'eac', 'cc'};
            for i = 1 : numel(eldes)
                for j = 1 : numel(cpts)
                    paramobj.(eldes{i}).(cpts{j}).thermalConductivity = tC;
                    paramobj.(eldes{i}).(cpts{j}).heatCapacity = hC;
                end
            end
            
            paramobj.SOC   = 0.5;
            paramobj.initT = 298.15;
            paramobj.J     = 0.001;
            paramobj.Ucut  = 2;
            
        end

    end
    
end
