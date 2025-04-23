
function sigma_g = sigma(omega, gamma,T,mu_c)

  e=1.6E-19;
  hbar=1.05E-34;
  kb=1.38E-23;

  sigma_intra = ((2*kb*T*e^2)*i)*(log(2*cosh(mu_c/(2*kb*T))))/(pi*hbar^2*(omega + i*gamma));
  sigma_inter = (e^2/(4*hbar))*(0.5 + atan((hbar*(omega + gamma*i) -2*mu_c)/(2*kb*T))/pi - i*log(((hbar*(omega+i*gamma)+2*mu_c)^2)/((hbar*(omega+i*gamma)-2*mu_c)^2 + (2*kb*T)^2))/(2*pi));
  sigma_g = sigma_inter + sigma_intra;

endfunction
