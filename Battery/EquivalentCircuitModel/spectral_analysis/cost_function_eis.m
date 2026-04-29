
function [total_error, Z_re_exp, Z_im_exp] = cost_function_eis(best_params)


Z_re_exp = plot_nyquist_ank.Z_real;
Z_im_exp = plot_nyquist_ank.Z_imag;


omega = logspace(-2, 4, 100)

Z_re_params = plot_nyquist_params.Z_real(best_params);
Z_im_params = plot_nyquist_params.Z_imag(best_params);

totalerror = sum((Z_re_exp(omega)-Z_re_params(omega)).^2) + sum((Z_im_exp(omega)-Z_im_params(omega))^2)