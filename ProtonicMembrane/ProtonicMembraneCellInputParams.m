classdef ProtonicMembraneCellInputParams < ComponentInputParams
    
    properties

        T
        
        Anode
        Cathode
        Electrolyte
        Control
        
        couplingTerms
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneCellInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
    
            pick = @(fd) pickField(jsonstruct, fd);
            paramobj.(an)    = ProtonicMembraneElectrodeInputParams(pick(an));
            paramobj.(ct)    = ProtonicMembraneElectrodeInputParams(pick(ct));
            paramobj.(elyte) = ProtonicMembraneElectrolyteInputParams(pick(elyte));
            paramobj.(ctrl)  = ProtonicMembraneElectrolyteInputParams(pick(ctrl));

            paramobj = mergeParameters(paramobj, {{'T'}       , ...
                                                  {elyte, 'T'}, ...
                                                  {an, 'T'}   , ...
                                                  {ct, 'T'}});
            
            paramobj.couplingTerms = {};
            
        end
        
    end
    
end
