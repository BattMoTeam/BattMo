function val = rampupControl(t, tup, inputI)

    rampupcase = 'sineup';        
    switch rampupcase
        
      case 'sineup'
        val = (t <= tup) .* sineup(0, inputI, 0, tup, t) + (t > tup) .* inputI;
        
      case 'linear'
        val = (t <= tup).*t./tup .* inputI +  (t > tup) .* inputI;
        
    end
    
end
