function [Z_real, Z_imag] = load_nyquist(params, omega)

    R_0 = params(1);
    R_1 = params(2);
    C_1 = params(3);
    R_2 = params(4);
    C_2 = params(5);
    
    
    
    Z_real = R_0 +    R_1 ./ (1+(R_1* C_1.*omega).^2)    +   R_2 ./ (1+(R_2* C_2.*omega).^2);
    Z_imag = (R_1*R_1*C_1.*omega ./ (1+(R_1* C_1.*omega).^2))  ...
               + (R_2*R_2*C_2.*omega ./ (1+(R_2* C_2.*omega).^2));

end