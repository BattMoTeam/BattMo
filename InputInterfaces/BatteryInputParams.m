classdef BatteryInputParams 
    
    properties
        
        % Global Grid
        G
        % Temperature and SOC
        T
        SOC
        % Input current
        J
        % Voltage cut
        Ucut
        
        %% parameters for the battery components
        % shortcut used here
        % ne : Negative electrode parameters (class ElectrodeInputParams)
        % pe : Positive electrode parameters (class ElectrodeInputParams)
        % elyte : Electrolyte (class ElectrolyteInputParams)
        ne;
        pe;
        elyte;        
                
        %% Coupling terms (describe the topological structure of the coupling)
        couplingTerms
        
    end
    
    methods
        
        function paramobj = BatteryInputParams()
            paramobj.ne = ElectrodeInputParams();
            paramobj.pe = ElectrodeInputParams();
            paramobj.elyte = ElectrolyteInputParams();
            parmobj.couplingTerms = {};
        end

        function paramobj = setup(paramobj, params)

            fdnames = {'ionName', ...
                       'ionFluxName', ...
                       'ionSourceName', ...
                       'ionMassConsName', ...
                       'ionAccumName'};
            
            paramobj = dispatchParams(paramobj, params, fdnames);
            
            fdnames = {'T', ...
                       'SOC', ...
                       'J', ...
                       'Ucut'};
            
            paramobj = dispatchParams(paramobj, params, fdnames);
                        
            fdnames = {'G', ...
                       'ne', ...
                       'pe', ...
                       'elyte', ...
                       'couplingTerms'};

            paramobj = dispatchParams(paramobj, params, fdnames);
            
        end


    end
    
end
