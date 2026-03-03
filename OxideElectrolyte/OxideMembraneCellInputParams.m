classdef OxideMembraneCellInputParams < ComponentInputParams
    
    properties

        T
        
        Anode
        Cathode
        Electrolyte
        Control
        
        couplingTerms

        dx
        farea
        
    end
    
    methods
        
        function inputparams = OxideMembraneCellInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
    
            pick = @(fd) pickField(jsonstruct, fd);
            inputparams.(an)    = OxideMembraneElectrodeInputParams(pick(an));
            inputparams.(ct)    = OxideMembraneElectrodeInputParams(pick(ct));
            inputparams.(elyte) = OxideMembraneElectrolyteInputParams(pick(elyte));
            inputparams.(ctrl)  = OxideMembraneControlInputParams(pick(ctrl));

            inputparams = mergeParameters(inputparams, {{'T'}       , ...
                                                  {elyte, 'T'}, ...
                                                  {an, 'T'}   , ...
                                                  {ct, 'T'}});
            
            inputparams = mergeParameters(inputparams, {{elyte, 'muEl0'}, ...
                                                  {ct, 'muEl0'}   , ...
                                                  {an, 'muEl0'}});

            inputparams = mergeParameters(inputparams, {{elyte, 'Keh'}, ...
                                                  {ct, 'Keh'}   , ...
                                                  {an, 'Keh'}});


            inputparams.couplingTerms = {};
            
        end
        
    end
    
end
