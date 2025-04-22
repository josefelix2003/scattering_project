%calculer sigma_g avec une fonction

function [omega, gamma1, gamma2] = omega_gamma_matrix (M, alpha, beta, d, k0, k1, k2, Z0, sigma_g, a)

  K=2*pi/d

  Z=zeros(2*M+1)
  M=k0*Z0*toeplitz_matrix(M, a, d, sigma_g)

  omega=[Z, M; M, Z];

  A=zeros(2*M+1)
  B=zeros(2*M+1)
  C=zeros(2*M+1)
  D=zeros(2*M+1)
  E=zeros(2*M+1)
  F=zeros(2*M+1)


  for p = 1:2*M+1
    n=p-(M+1) #p matrix index and n true channel index
    alpha_n=alpha+n*K
    gamma_1n=sqrt(k1^2 - alpha_n^2 - beta^2)
    gamma_2n=sqrt(k2^2 - alpha_n^2 - beta^2)

    A(p, p) = alpha_n*beta/gamma_1n
    B(p, p) = gamma_1n + (beta^2/gamma_1n)
    C(p, p) = gamma_1n + (alpha_n^2/gamma_1n)

    D(p, p) = alpha_n*beta/gamma_2n
    E(p, p) = gamma_2n + (beta^2/gamma_2n)
    F(p, p) = gamma_2n + (alpha_n^2/gamma_2n)



  endfor

  gamma1=[A, B; C, A]
  gamma2=[D, E; F, D]


endfunction
