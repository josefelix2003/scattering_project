
function T = toeplitz_matrix (M, a, d, sigma_g)
  c_vect=[];
  r_vect=[];
  for n = 0:2*M
    c_vect(end+1) = get_sigma_coeff(n, a, d, sigma_g);
    r_vect(end+1) = get_sigma_coeff(-n, a, d, sigma_g);

  endfor

  T = toeplitz(c_vect,r_vect)

endfunction
