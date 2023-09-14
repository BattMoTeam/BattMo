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
        
        function paramobj = OxideMembraneCellInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
            an    = 'Anode';
            ct    = 'Cathode';
            elyte = 'Electrolyte';
            ctrl  = 'Control';
    
            pick = @(fd) pickField(jsonstruct, fd);
            paramobj.(an)    = OxideMembraneElectrodeInputParams(pick(an));
            paramobj.(ct)    = OxideMembraneElectrodeInputParams(pick(ct));
            paramobj.(elyte) = OxideMembraneElectrolyteInputParams(pick(elyte));
            paramobj.(ctrl)  = OxideMembraneControlInputParams(pick(ctrl));

            paramobj = mergeParameters(paramobj, {{'T'}       , ...
                                                  {elyte, 'T'}, ...
                                                  {an, 'T'}   , ...
                                                  {ct, 'T'}});
            
            paramobj = mergeParameters(paramobj, {{ct, 'muEl0'}, ...
                                                  {an, 'muEl0'}});

            paramobj = mergeParameters(paramobj, {{elyte, 'Keh'}, ...
                                                  {ct, 'Keh'}   , ...
                                                  {an, 'Keh'}});


            paramobj.couplingTerms = {};
            
        end
        
    end
    
end
