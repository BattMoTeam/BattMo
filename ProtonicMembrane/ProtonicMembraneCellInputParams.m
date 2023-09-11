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
            paramobj.(an)    = ProtonicMembraneAnodeInputParams(pick(an));
            paramobj.(ct)    = ProtonicMembraneCathodeInputParams(pick(ct));
            paramobj.(elyte) = ProtonicMembraneElectrolyteInputParams(pick(elyte));
            paramobj.(ctrl)  = ProtonicMembraneControlInputParams(pick(ctrl));

            paramobj = mergeParameters(paramobj, {{'T'}       , ...
                                                  {elyte, 'T'}, ...
                                                  {an, 'T'}   , ...
                                                  {ct, 'T'}});
            
            paramobj = mergeParameters(paramobj, {{elyte, 'EO2_0'}, ...
                                                  {an, 'Eocp'}});

            paramobj = mergeParameters(paramobj, {{elyte, 'EH2_0'}, ...
                                                  {ct, 'Eocp'}});

            paramobj = mergeParameters(paramobj, {{elyte, 'SU'}, ...
                                                  {an, 'SU'}});

            paramobj = mergeParameters(paramobj, {{elyte, 'pH2O_in'}, ...
                                                  {an, 'pH2O_in'}});

            paramobj.couplingTerms = {};
            
        end
        
    end
    
end
