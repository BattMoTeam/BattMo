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
        
        function inputparams = ProtonicMembraneCellInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
    
            pick = @(fd) pickField(jsonstruct, fd);
            inputparams.(an)    = ProtonicMembraneAnodeInputParams(pick(an));
            inputparams.(ct)    = ProtonicMembraneCathodeInputParams(pick(ct));
            inputparams.(elyte) = ProtonicMembraneElectrolyteInputParams(pick(elyte));
            inputparams.(ctrl)  = ProtonicMembraneControlInputParams(pick(ctrl));

            inputparams = mergeParameters(inputparams, {{'T'}       , ...
                                                  {elyte, 'T'}, ...
                                                  {an, 'T'}   , ...
                                                  {ct, 'T'}});
            
            inputparams = mergeParameters(inputparams, {{elyte, 'Ptot'}, ...
                                                  {an, 'Ptot'}   , ...
                                                  {ct, 'Ptot'}});

            inputparams = mergeParameters(inputparams, {{elyte, 'SU'}, ...
                                                  {an, 'SU'}});

            inputparams = mergeParameters(inputparams, {{elyte, 'steam_ratio'}, ...
                                                  {an, 'steam_ratio'}});

            inputparams = mergeParameters(inputparams, {{elyte, 'E_0'}, ...
                                                  {an, 'E_0'}});
            
            inputparams.couplingTerms = {};
            
        end
        
    end
    
end