classdef Electrolyte < ElectroChemicalComponent
    
    properties
        
        sp
        compnames
        ncomp
        
        volumeFraction
        
    end

    methods
        
        function model = Electrolyte(paramobj)
            
            model = model@ElectroChemicalComponent(paramobj);
            
            fdnames = {'sp', ...
                       'compnames', ...
                       'ncomp'};
            model = dispatchParams(model, paramobj, fdnames);
            
            
        end

        function state = updateChemicalCurrent(model, state)
            error('virtual function')
        end
        
    end
end

