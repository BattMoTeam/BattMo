function [I, alpha, beta] = controlfunc(time, Imax, tswitch, T, varargin)

    opt = struct('order', 'alpha-first');
    opt = merge_options(opt, varargin{:});
    
    switch opt.order
      case 'alpha-first'
        if time <= tswitch
            I     = 0;
            alpha = time/tswitch;
            beta  = 0;
        else
            I     = (time - tswitch)/(T - tswitch)*Imax;
            alpha = 1;
            beta  = (time - tswitch)/(T - tswitch);
        end
      case 'alpha-first-beta-constant'
        if time <= tswitch
            I     = 0;
            alpha = time/tswitch;
            beta  = 0;
        else
            I     = (time - tswitch)/(T - tswitch)*Imax;
            alpha = 1;
            beta  = 1;
        end
      case 'alpha-equal-beta'
        if time <= tswitch
            I     = time/tswitch*Imax;
            alpha = 0;
            beta  = 0;
        else
            I     = Imax;
            alpha = (time - tswitch)/(T - tswitch);
            beta  = alpha;
        end        
      case 'I-first'
        if time <= tswitch
            I     = time/tswitch*Imax;
            alpha = 0;
            beta  = time/tswitch;
        else
            I     = Imax;
            alpha = (time - tswitch)/(T - tswitch);
            beta = 1;
        end
      case 'alpha-only'
        % tswitch is not used. This case is used in OxideElectrolyte where there is no Butler-Volmer at boundary.
        I = Imax;
        alpha = time/T;
      otherwise
        error('order not recognized');
    end
    
end
