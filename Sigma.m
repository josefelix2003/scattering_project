function Sigma = Sigma(w, damp, T, mu_c)
  % Calcul rapide de la conductivité complexe (Sigma)
  % w : fréquence angulaire en rad/s
  % damp : terme d’amortissement
  % T : température (K)
  % mu_c : potentiel chimique (J)

  % Constantes utiles
  e = 1.602e-19;
  hbar = 1.054e-34;
  kb = 1.381e-23;

  % Valeurs par défaut (si rien passé)
  if nargin < 4
    mu_c = 0.03 * e; % ~30 meV
  end
  if nargin < 3
    T = 290; % un peu en dessous de 300
  end
  if nargin < 2
    damp = 0.08; % légèrement plus faible
  end
  if nargin < 1
    w = logspace(11, 14.5, 250); % plage raisonnable
  end

  % Partie intrabande (valable surtout pour basses fréquences)
  num_intra = 2 * kb * T * e^2 * 1i * log(2 * cosh(mu_c / (2 * kb * T)));
  denom_intra = pi * hbar^2 * (w + 1i * damp);
  sigma_intra = num_intra ./ denom_intra;

  % Partie interbande (effets à haute fréquence)
  omega_c = hbar * (w + 1i * damp);
  x = (omega_c - 2 * mu_c) ./ (2 * kb * T);
  term1 = 0.5 + atan(x) / pi;
  term2 = log(((omega_c + 2*mu_c).^2) ./ ((omega_c - 2*mu_c).^2 + (2*kb*T)^2));
  sigma_inter = (e^2 / (4 * hbar)) * (term1 - 1i * term2 / (2 * pi));

  % Somme totale
  Sigma = sigma_intra + sigma_inter;

  % Tracés : 2 styles pour mieux voir
  figure;
  loglog(w, abs(Sigma), 'b', 'LineWidth', 1.3);
  grid on;
  xlabel('\omega (rad/s)');
  ylabel('|Σ(ω)|');
  title('Conductivité (log-log)');

  figure;
  plot(w, abs(Sigma), 'k', 'LineWidth', 1.2);
  grid on;
  xlabel('\omega (rad/s)');
  ylabel('|Σ(ω)|');
  title('Conductivité (linéaire)');
end

% Lancer avec les valeurs par défaut
Sigma();


