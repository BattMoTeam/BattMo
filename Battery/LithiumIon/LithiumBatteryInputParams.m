classdef LithiumBatteryInputParams < BatteryInputParams
% reference
% @misc{Latz_2016,
% 	doi = {10.1515/nano.bjneah.6.102},
% 	url = {https://doi.org/10.1515%2Fnano.bjneah.6.102},
% 	year = 2016,
% 	month = {jul},
% 	publisher = {De Gruyter},
% 	author = {Arnulf Latz and Jochen Zausch},
% 	title = {Multiscale modeling of lithium ion batteries: thermal aspects}
% }
    
    methods
        
        function paramobj = LithiumBatteryInputParams()
            
            paramobj = paramobj@BatteryInputParams();

            paramobj.ne = LithiumElectrodeInputParams();
            paramobj.pe = LithiumElectrodeInputParams();

            % We set interdiffusion coefficient here
            paramobj.ne.eac.InterDiffusionCoefficient = 1e-10;
            paramobj.pe.eac.InterDiffusionCoefficient = 1e-10;
            
            % Once the electrode are instantiated, we fill in the active material inputs
            paramobj.ne.eac.am = GraphiteInputParams();
            paramobj.pe.eac.am = NMC111InputParams();
            
            paramobj.elyte = orgLiPF6InputParams();
            
            % Values from Latz (reference Latz_2016)
            tC = 0.006;  % Thermal Conductivity (in W/(cm K))
            rho = 0.001; % Density (in kg/(cm^3))
            hC = 4180;   % Heat Capacity (in J/(kg K))
            
            % volumetric heat capacity (in SI unit: J/(m^3 K))
            hC = rho*hC*1e6;
            % thermal conducivity (in SI unit: W/(m K))
            tC = tC*1e-2;
            
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
