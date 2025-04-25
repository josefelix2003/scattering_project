function [I, J, R, T] = IJRT_vectors (Ixy, Jxy, S, M, d, alpha, beta, k1, k2)
  K=2*pi/d;

  I_J = [Ixy;Jxy];
  R_T = S*I_J;

  R=R_T(1:4*M+2);
  T=R_T(4*M+3:8*M+4);

  I = Ixy;
  J = Jxy;

  for px = 1:2*M+1 %px: x vecteur index
    py = px+1+2*M;  %py: y vecteur index
    n = px-1-M;

    alpha_n=alpha+n*K;
    gamma_1n=sqrt(k1^2 - alpha_n^2 - beta^2);
    gamma_2n=sqrt(k2^2 - alpha_n^2 - beta^2);


    I(end+1)=-(alpha_n*I(px) + beta*I(py))/gamma_1n;
    J(end+1)=(alpha_n*J(px) + beta*J(py))/gamma_2n;
    R(end+1)=(alpha_n*R(px) + beta*R(py))/gamma_1n;
    T(end+1)=-(alpha_n*T(px) + beta*T(py))/gamma_2n;


  endfor


endfunction
