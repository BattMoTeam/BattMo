function [I, alpha, beta] = pmControlFunc(time, targetI, tswitch, T, varargin)
    
    opt = struct('order', 'alpha-first');
    opt = merge_options(opt, varargin{:});
    
    switch opt.order
      case 'alpha-first'
        if time <= tswitch
            I     = 0;
            alpha = time/tswitch;
            beta  = 0;
        else
            I     = (time - tswitch)/(T - tswitch)*targetI;
            alpha = 1;
            beta  = (time - tswitch)/(T - tswitch);
        end
      case 'alpha-first-beta-constant'
        if time <= tswitch
            I     = 0;
            alpha = time/tswitch;
            beta  = 0;
        else
            I     = (time - tswitch)/(T - tswitch)*targetI;
            alpha = 1;
            beta  = 1;
        end
      case 'alpha-equal-beta'
        if time <= tswitch
            I     = time/tswitch*targetI;
            alpha = 0;
            beta  = 0;
        else
            I     = targetI;
            alpha = (time - tswitch)/(T - tswitch);
            beta  = alpha;
        end        
      case 'I-first'
        if time <= tswitch
            I     = time/tswitch*targetI;
            alpha = 0;
            beta  = time/tswitch;
        else
            I     = targetI;
            alpha = (time - tswitch)/(T - tswitch);
            beta = 1;
        end
      case 'alpha-only'
        % tswitch is not used. This case is used in OxideElectrolyte where there is no Butler-Volmer at boundary.
        I = targetI;
        alpha = time/T;
      otherwise
        error('order not recognized');
    end
    
end
