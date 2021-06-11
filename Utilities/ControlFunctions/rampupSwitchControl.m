function [val, ctrltype] = rampupSwitchControl(t, tup, I, E, inputI, inputE)
% CURRENTSOURCE Summary of this function goes here
%   Detailed explanation goes here

    if (E > inputE) | (t < tup)
        
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
