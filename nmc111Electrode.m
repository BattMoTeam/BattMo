classdef nmc111Electrode < FvModel
    %NMC111ELECTRODE Summary of this class goes here
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
        cei     % Cathode-electrolyte interphase (CEI) object
        elyte   % Liquid electrolyte data structure
        
        sigmaeff
        
        % State properties
        E       % Electric potential,   [V]
        eta     % Overpotential,        [V]
        j       % Current density,      [A/m2]
        R
        T       % Temperature
        
        % Mesh properties
        X
        Xb
        N
        dombin
        
    end
    
    methods
        function obj = nmc111Electrode(compname, SOC, T, G, cells)
            
            obj = obj@FvModel(compname);
            
            obj.am = nmc111AM(SOC, T);
            obj.bin = ptfe();
            
            obj.am.eps = 0.8;
            
            obj.eps = obj.am.eps + obj.bin.eps;
                    
            obj.t = 10e-6;
            obj.T = T;
            
            obj.thermodynamics()
            obj.E = obj.am.OCP;
        
            obj.Grid = genSubGrid(G, cells);
            
            obj.varnames = {'phi', 'Li'};
            nc = obj.Grid.cells.num;
            obj.varsizes = [nc, nc];
            
        end
        
        function thermodynamics(obj)
            %THERMODYNAMICS Summary of this method goes here
            %   Detailed explanation goes here
           
            
            % Electrochemical Thermodynamics
            obj.am.equilibrium();
            
            
        end
        
        function reactBV(obj, phiElyte)
            
            obj.thermodynamics();
            
            obj.eta    =   -(obj.am.phi - ...
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

