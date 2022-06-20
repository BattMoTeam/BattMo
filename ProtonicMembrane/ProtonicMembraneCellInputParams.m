classdef ProtonicMembraneCellInputParams < ComponentInputParams
    
    properties
        
        Anode
        Cathode
        Electrolyte
        couplingTerms
    end
    
    methods
        
        function paramobj = ProtonicMembraneCellInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            
            pick = @(fd) pickField(jsonstruct, fd);
            paramobj.(an)    = ProtonicMembraneElectrodeInputParams(pick(an));
            paramobj.(ct)    = ProtonicMembraneElectrodeInputParams(pick(ct));
            paramobj.(elyte) = ProtonicMembraneElectrolyteInputParams(pick(elyte));
            
            paramobj.couplingTerms = {};
        end
        
    end
    
end
