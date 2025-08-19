function [OCP, dUdT] = computeOCP_LFP_Gerver2011(c, T, cmax)

    error("function is not compatible with new function interface");
    
    OCP = 3.41285712e+00 ...
          - 1.49721852e-02 * c/cmax ...
          + 3.54866018e+14 * exp(-3.95729493e+02 * c/cmax) ...
          - 1.45998465e+00 * exp(-1.10108622e+02 * (1 - c/cmax));
    
    dUdT = 0;
    
end
