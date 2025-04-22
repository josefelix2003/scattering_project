
function sigma_p = get_sigma_coeff (p, a, d, sigma_g)
  if (p==0)
    sigma_p = sigma_g * a / d;
  else
    K = 2*pi/d;
    sigma_p = (exp(-i*K*p*a) - 1)*( i*sigma_g/(2*p*pi) );
  endif

endfunction
