function sigma_g = sigma_g(w, damp, T, mu_c)
  % sigma_g = sigma_g(w, damp, T, mu_c)
  % Calcule la conductivité totale sigma_g en fonction de la pulsation temporelle w,
  % de l'amortissement damp, de la température T et du potentiel chimique mu_c.
  %
  % w : la pulsation temporelle (rad/s)
  % damp : l'amortissement (sans dimension)
  % T : la température (en Kelvin)
  % mu_c : le potentiel chimique (en Joules)

  e = 1.602e-19; % charge élémentaire (Coulomb)
  hbar = 1.054e-34; % constante de Planck réduite (J.s)
  kb = 1.381e-23; % constante de Boltzmann (J/K)

  % Paramètres
  damp = 0.1;        % Amortissement (sans dimension)
  T = 300;           % Température en Kelvin
  mu_c = 0.026 * 1.602e-19; % Potentiel chimique en Joules (environ 0.026 eV)

% Générer une plage de pulsations temporelles (omega)
  w = logspace(10, 15, 200); % Balayage logarithmique de 10^10 à 10^15 rad/s

  sigma_intra = ((2*kb*T*e^2)*1i)*(log(2*cosh(mu_c/(2*kb*T)))) ./ (pi*hbar^2*(w + 1i*damp));
  % Conductivité intrabande dominante à basse fréquence.
  sigma_inter = (e^2/(4*hbar))*(0.5 + atan((hbar*(w + damp*1i) -2*mu_c)/(2*kb*T))/pi - 1i*log(((hbar*(w+1i*damp)+2*mu_c).^2)./((hbar*(w+1i*damp)-2*mu_c).^2 + (2*kb*T)^2))/(2*pi));
  % Conductivité interbande, dominante à haute fréquence.
  sigma_g = sigma_inter + sigma_intra;
  % Conductivité totale.

  % Tracer le module de sigma_g en fonction de omega (échelle logarithmique)
  figure;
  loglog(w, abs(sigma_g), 'b-', 'linewidth', 1.5);
  grid on;
  xlabel('\omega (rad/s)');
  ylabel('|\sigma_g(\omega)|');
  title('Module de la conductivité en fonction de la pulsation temporelle (échelle logarithmique)');



  % Tracer le module de sigma_g en fonction de omega (échelle lineaire)
  figure;
  plot(w,abs(sigma_g))
  grid on;
  xlabel('\omega (rad/s)');
  ylabel('|\sigma_g(\omega)|');
  title('Module de la conductivité en fonction de la pulsation temporelle (échelle logarithmique)');
endfunction



% Calculer et afficher la conductivité
sigma_g(w, damp, T, mu_c);
