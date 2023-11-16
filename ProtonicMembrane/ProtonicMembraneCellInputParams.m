classdef ProtonicMembraneCellInputParams < ComponentInputParams
    
    properties

        T
        
        Anode
        Cathode
        Electrolyte
        Control
        
        couplingTerms

        dx
        faceArea
        
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
            
            paramobj = mergeParameters(paramobj, {{elyte, 'Ptot'}, ...
                                                  {an, 'Ptot'}   , ...
                                                  {ct, 'Ptot'}});

            paramobj = mergeParameters(paramobj, {{elyte, 'SU'}, ...
                                                  {an, 'SU'}});

            paramobj = mergeParameters(paramobj, {{elyte, 'steam_ratio'}, ...
                                                  {an, 'steam_ratio'}});

            paramobj = mergeParameters(paramobj, {{elyte, 'E_0'}, ...
                                                  {an, 'E_0'}});
            
            paramobj.couplingTerms = {};
            
        end
        
    end
    
end
