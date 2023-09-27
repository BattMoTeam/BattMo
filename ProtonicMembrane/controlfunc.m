function [I, alpha] = controlfunc(time, Imax, tswitch, T, varargin)

    opt = struct('order', 'alpha-first');
    opt = merge_options(opt, varargin{:});
    
    switch opt.order
      case 'alpha-first'
        if time <= tswitch
            I = 0;
            alpha = time/tswitch;
        else
            I = (time - tswitch)/(T - tswitch)*Imax;
            alpha = 1;
        end
      case 'I-first'
        if time <= tswitch
            I     = time/tswitch*Imax;
            alpha = 0;
        else
            I     = Imax;
            % alpha = 0;
            alpha = (time - tswitch)/(T - tswitch);
        end
      case 'alpha-only'
        % tswitch is not used. This case is used in OxideElectrolyte where there is no Butler-Volmer at boundary.
        I = Imax;
        alpha = time/T;
      otherwise
        error('order not recognized');
    end
    
end
