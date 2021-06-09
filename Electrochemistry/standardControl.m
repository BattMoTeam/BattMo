function [val, ctrltype] = standardControl(t, tup, tlim, I, E, inputI, inputE)
% CURRENTSOURCE Summary of this function goes here
%   Detailed explanation goes here

    if E > inputE
        
        % We control with current
        ctrltype = 'I';
        
        rampupcase = 'sineup';        
        switch rampupcase
            
          case 'sineup'
            val = (t <= tup) .* sineup(0, inputI, 0, tup, t) + (t > tup && t <= tlim) .* inputI;
            
          case 'linear'
            val = (t <= tup).*t./tup .* inputI +  (t > tup && t <= tlim) .* inputI;
            
        end
    
    else
        
        % We control with current
        ctrltype = 'E';

        val = inputE;
        
    end
    
end
