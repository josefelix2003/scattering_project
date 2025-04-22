
function S = scatttering_matrix (gamma1, gamma2, omega, M)
  A=eye(4*M+2)
  B=-A

  C=gamma1
  D=gamma2+omega
  E=gamma2-omega

  M1=[A, B; C, D]
  M2=[B, A; C, E]

  S=inv(M1)*M2

endfunction
