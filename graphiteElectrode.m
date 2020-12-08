classdef graphiteElectrode < FvModel
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = physicalConstants();
        
        % Design properties
        t       % Thickness,        [m]
        eps     % Volume fraction,  [-]
        void    % Porosity,         [-]
        
        % Material properties
        am      % Active material object
        bin     % Binder object
        ca      % Conducting additive object
        sei     % Solid-electrolyte interphase (SEI) object
        elyte   % Liquid electrolyte data structure
        
        % Effective conductivity
        sigmaeff
        
        % State properties
        E       % Electric potential,   [V]
        eta     % Overpotential,        [V]
        j       % Current density,      [A/m2]
        R       % Reaction Rate
        T       % Temperature
        
        % Mesh properties (updated before simulation)
        X
        Xb
        N
        dombin
                
    end
    
    methods
        function obj = graphiteElectrode(SOC, T, dims, sizes)
            %UNTITLED6 Construct an instance of this class
            %   Detailed explanation goes here
            obj.am  = graphiteAM(SOC, T);
            obj.bin = ptfe();
            obj.sei = seiAM();
            
            obj.am.eps = 0.8;
            
            obj.eps =   obj.am.eps ...
                        + obj.bin.eps ...
                        + obj.sei.eps;
            
            obj.t = 10e-6;
            
            obj.T = T;
            
            obj.thermodynamics()
            obj.E = obj.am.OCP;
            
            obj.Grid = cartGrid(dims, sizes);
            
        end
        
        function thermodynamics(obj)
            %THERMODYNAMICS Summary of this method goes here
            %   Detailed explanation goes here
           
            
            % Electrochemical Thermodynamics
            obj.am.equilibrium();
            
            
        end
        
        function reactBV(obj, phiElyte)
            
            obj.thermodynamics();
            
            obj.eta    =   (obj.am.phi - ...
                                phiElyte - ...
                                obj.am.OCP);
                                    
            obj.R   =   obj.am.Asp .* butlerVolmer(   obj.am.k .* 1 .* obj.con.F, ...
                                        0.5, ...
                                        1, ...
                                        obj.eta, ...
                                        obj.T ) ./ (1 .* obj.con.F);
                                    
        end
        
    end
end

