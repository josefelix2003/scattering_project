hbar = 1.05^-34;
k_b = 1.38^-23;
e = 1.602^-19;
T = 300;
mu_c = 1;
omega = 1e13;
gamma = 1e12;

function [sigma_intra] = sigma(omega, gamma, e, hbar, k_b, T, mu_c)
  sigma_intra = (2*i * e*e * k_b * T) / (pi * hbar^2 * (omega + 1*i * gamma)) * log(2 * cosh(mu_c / (2 * k_b * T)));
endfunction

sigma_intra = sigma(omega, gamma, e, hbar, k_b, T, mu_c);

disp(sigma_intra)
