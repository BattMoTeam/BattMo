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

        function paramobj = setup(paramobj, params)
        % params struct should contain valid fields for ComponentInputParams,
        %
        % valid fields for the methods (see implementation of those methods)
        %
        %
        % and fields
        %
        % - T
        % - SOC
        % - J
        % - Ucut
        %
        % - G
        % - ne
        % - pe            
        % - elyte
        % - couplingTerms
            
            paramobj = paramobj.setupVariousParams(params);
            paramobj.G  = getparam(params, G);
            paramobj.ne = getparam(params, 'ne');
            paramobj.pe = params.pe;
            paramobj.elyte = params.elyte;
            paramobj.coupltingTerms = params.couplingTerms;
        end

        function paramobj = setupVariousParams(paramobj, params)
            paramobj.T    = params.T;
            paramobj.SOC  = params.SOC;
            paramobj.J    = params.J;
            paramobj.Ucut = params.Ucut;
        end


    end
    
end
