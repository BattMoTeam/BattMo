function [val, ctrltype] = rampupSwitchControl(t, tup, I, E, inputI, inputE)
% CURRENTSOURCE Summary of this function goes here
%   Detailed explanation goes here

    if(inputI<0)
        %notuseVolt = ((E - inputE)<1e-3);
        %notuseVolt = not(useVolt);
        useCurrent = (E < inputE);
    else
        %notuseVolt = ((E - inputE)>-1e-3);
        useCurrent = (E > inputE);
        %useVolt = not(useVolt);
    end
    %if not(useVolt) || (t < tup)
    if  useCurrent | (t < tup)   
        % We control with current
        ctrltype = 'I';
        
        rampupcase = 'sineup';        
        switch rampupcase
            
          case 'sineup'
            val = (t <= tup) .* sineup(0, inputI, 0, tup, t) + (t > tup) .* inputI;
            
          case 'linear'
            val = (t <= tup).*t./tup .* inputI +  (t > tup ) .* inputI;
            
        end
    
    else
        
        % We control with current
        ctrltype = 'E';

        val = inputE;
        
    end
    
end
