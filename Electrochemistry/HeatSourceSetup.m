classdef HeatSourceSetup

    properties

        sourceTerms
        times

    end

    methods

        function hss = HeatSourceSetup(sourceTerms, times)

            hss.sourceTerms = sourceTerms;
            hss.times       = times;
            
        end

        function sourceTerm = eval(hss, time)

            times       = hss.times;
            sourceTerms = hss.sourceTerms;
            
            N = numel(times);
            
            ind = interp1(times, (1 : N)', time);

            ind0 = floor(ind);
            ind1 = ind0 + 1;

            if ind1 > N
                sourceTerm = sourceTerms{ind0};
            else
                sourceTerm = (ind1 - ind)*sourceTerms{ind0} + (ind - ind0)*sourceTerms{ind1};
            end

        end
        
    end
    
    
end

