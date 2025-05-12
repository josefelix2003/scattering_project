function results = run_simulation ()
  clear
  config %On charge le fichier config.m pour parametrer les variables de la simulation.

function results = compute_values (in, out) % "in" les input et "out" les output extraits de config.

  #--------------------Constants-----------------------
  c=3E8; %Vitesse de la lumière dans le vide (m/s).
  hbar=1.05E-34; %Constante de Planck réduite (J·s).
  permittivity0=8.8541E-12; %Permittivité du vide (F/m). F : farad
  permeability0=1.2566E-6; %Perméabilité du vide (H/m) ou (T*m/A).
  Z0=sqrt(permeability0/permittivity0); %Impédance du milieu (Ω).
  damping_factor=0.658E-22/hbar; %Facteur d'amortissement (1/s).
  kb=1.38E-23; % constante de Boltzmann (J/K).
  e=1.6E-19; %charge élémentaire (Coulomb)


  #----------------------------------------------------
  #-----------------Derived quantities-----------------
  k0=2*pi/in.lambda; %Vecteur d'onde dans le vide (rad/m).
  n=sqrt(in.permittivityr1);
  w=2*pi*c/(n*in.lambda); %Pulsation (rad/s).
  sigma_g=sigma(w, damping_factor, in.T, in.mu_c); %Conductivité du graphène (S/m).
  k1=k0*sqrt(in.permittivityr1*in.permeabilityr1); %Vecteur d'onde milieu 1 (rad/m).
  k2=k0*sqrt(in.permittivityr2*in.permeabilityr2); %Vecteur d'onde milieu 2 (rad/m).

                  %Incident wave vector
  alpha=k1*sin(in.theta)*cos(in.phi); %Composante x du vecteur k (rad/m).
  beta=k1*sin(in.theta)*sin(in.phi);  %Composante y vecteur k (rad/m).
  gamma=k1*cos(in.theta);             %Composante z vecteur k (rad/m).
  #----------------------------------------------------
  #-------------Computing Scattering Matrix------------
              %Omega, Gamma1, Gamma2 Matrices
  [omega, gamma1, gamma2] = omega_gamma_matrix(in.M, alpha, beta, in.d, k0, k1, k2, Z0, sigma_g, in.a); %Appel de la fonction qui calcule les matrices omega, gamma1, gamma2.
  % M : le nombre de strips de graphène. % alpha : Composante x vecteur k (rad/m).  % beta : Composante y vecteur k (rad/m). % d : Période du réseau (m)  % sigma_g : La conductivité sigma du graphène (S/m) % a : Largeur d'une bande de graphène (m).
  % : Pour le reste cf. Derived quantities
                    %Scattering Matrix
  S=compute_scattering_matrix(gamma1, gamma2, omega, in.M); %Appel de la fonction qui calcule la matrice S.
  % omega, gamma1, gamma2 sont les matrices qui permettent le calcule de S et calculées avec la fonction omega_gamma_matrix.

  #----------------------------------------------------
  #-----------------Compute R, T vectors---------------
  [I, J, R, T]=IJRT_vectors(in.Ixy, in.Jxy, S, in.M, in.d, alpha, beta, k1, k2); %Appel de la fonction qui calcule les 3 composants des vecteurs I, J, R, T.
  %Ixy et Jxy les deux premières composantes de I et J (x et y) ; calculés à partir de la matrice S à l'ordre de troncature M. ; chacun des vecteurs sous la forme [Ix, Iy, Iz] et de taille (2*M+1)*3.
    % : Pour le reste cf. Derived quantities et Incident wave vector
  #----------------------------------------------------
  #--------------------Absorption----------------------
  A=compute_absorption(R, T, in.M, alpha, beta, in.d, k1, k2); %Apelle de la fonction qui calcule l'absorption lors de l'interaction avec le graphène.
  % R et T contiennent trois composantes (selon x, y, z).
  #----------------------------------------------------
  #-------------Building Results structure-------------
  results=struct(); % Création d'une structure vide pour stocker les résultats

  if out.absorption
    results.absorption = A; % Stocke l'absorption si demandé
  endif
  if out.I
    results.I = I; % Stocke le vecteur I si demandé
  endif
  if out.J
    results.J = J; % Stocke le vecteur J si demandé
  endif
  if out.R
    results.R = R; % Stocke le vecteur R si demandé
  endif
  if out.T
    results.T = T; % Stocke le vecteur T si demandé
  endif
  if out.scattering_matrix
    results.scattering_matrix = S; % Stocke la matrice S si demandé
  endif
  if out.conductivity
    results.conductivity = sigma_g; % Stocke la conductivité si demandé
  endif

endfunction %Retourne la structure Results.

%Balayage (sweep) si l'option sweep.enable est vraie dans le fichier config :

  if sweep.enable
    sweep_var=sweep.input_variable; %Paramètre entrant à balayer.
    output_var=sweep.output_variable; %Paramètre sortant à balayer.
    min_val=sweep.min; %%Valeur minimale du paramètre.
    max_val=sweep.max; %Valeur maximale du paramètre.
    step_val=sweep.step; %pas entre deux valeurs de la variable de balayage.
    xaxis=zeros(0,1); %Création des axes (utilisés pour faire le plot) pour les input .
    yaxis=zeros(0,1); %Création des axes (utilisés pour faire le plot) pour les output.
    counter=0; %Compteur d'itération.

    for current_val = min_val:step_val:max_val
      counter+=1; % Incrémentation du compteur d'itérations.
      in.(sweep_var)=current_val; % Mise à jour de la valeur de la variable de balayage.
      sim_results = compute_values(in, out);
      xaxis(counter,1)= current_val;  % Enregistrement de la valeur du paramètre balayé dans l'axe des x.
      yaxis(counter,1)= sim_results.(output_var);  % Enregistrement de la valeur de sortie dans l'axe des y.
    endfor
    results=struct();
    results.xaxis=xaxis; % Axe des x : paramètre balayé.
    results.yaxis=yaxis; % Axe des y : paramètre de sortie.

  #----------------------------------------------------
  #---------------Building Results structure-----------

   % Obtient le chemin absolu du fichier actuel (script ou fonction).
    current_file_path = mfilename('fullpath');
    [folder_path, ~, ~] = fileparts(current_file_path);

    % Crée le chemin vers le dossier "results" situé dans le même dossier que le script.
    result_folder = fullfile(folder_path, 'results');

    % Crée le dossier "results" s’il n’existe pas déjà.
    if ~exist(result_folder, 'dir')
        mkdir(result_folder);
    end

    % Spécifie le chemin complet pour enregistrer les résultats dans un fichier texte.
    fichier = fullfile(result_folder, 'results.txt');
    fid = fopen(fichier, 'w');

    % Vérifie que le fichier est bien ouvert.
    if fid == -1
        error('Impossible d''ouvrir le fichier.');
    end

    % Écrit l'en-tête dans le fichier texte.
    fprintf(fid, 'Simulation Results:\n\n');
    fprintf(fid, '%s\t%s\n', sweep_var, output_var);

    % Écrit les résultats ligne par ligne.
    for i = 1:length(xaxis)
        fprintf(fid, '%.4e\t%.4e\n', xaxis(i), yaxis(i));
    end

    % Ferme le fichier
    fclose(fid);

    % Si la condition pour afficher le plot dans le fichier config est vraie.
    if sweep.plot
        % Crée le plot
        plot(xaxis, yaxis);
        xlabel(sweep_var, 'FontSize', 16, 'FontName', 'Times');
        ylabel(output_var, 'FontSize', 16, 'FontName', 'Times');
        grid on;
        set(gca, 'FontSize', 12, 'FontName', 'Times');
        axis tight;
        box on;

        % Génère le titre en fonction des variables.
        title_str = sprintf('%s vs. %s', output_var, sweep_var);
        title(title_str, 'FontSize', 16, 'FontWeight', 'bold');

        % Définit le chemin du fichier pour enregistrer la figure.
        plot_filename = fullfile(result_folder, 'plot.png');

        % Sauvegarde la figure en PNG.
        print(gcf, plot_filename, '-dpng');
    end


    return

  #----------------------------------------------------

  endif

  if ~sweep.enable % Si le balayage est désactivé (calcule simple, pas de variation de paramètre).

  results=compute_values(in, out); % Calcul des résultats en appelant la fonction compute_values.

  current_file_path = mfilename('fullpath'); % Récupère le chemin complet du fichier actuel.
  [folder_path, ~, ~] = fileparts(current_file_path);

    % Crée le chemin vers le dossier "resultats".
  result_folder = fullfile(folder_path, 'results');

    % Crée le dossier s’il n’existe pas.
  if ~exist(result_folder, 'dir')
    mkdir(result_folder);
  end

  fichier = fullfile(result_folder, 'results.txt');   % Choix du nom du fichier de sortie.
  fid = fopen(fichier, 'w'); % Ouverture du fichier en écriture.

  if fid == -1
    error('Impossible d''ouvrir le fichier.'); % Si le fichier ne peut pas être ouvert.
  end

  fprintf(fid, 'Simulation Results:\n\n');   % Nom du header du fichier.
  % Sauvegarde des résultats selon le choix des sorties.
  if out.absorption
    fprintf(fid, 'Absorption = %.4e\n', results.absorption);
  endif
  if out.I
    save_vector(fid, results.I, "I"); % Sauvegarde du vecteur I.
    fprintf(fid, '\n');
  endif
  if out.J
    save_vector(fid, results.J, "J"); % Sauvegarde du vecteur J.
    fprintf(fid, '\n');
  endif
  if out.R
    save_vector(fid, results.R, "R"); % Sauvegarde du vecteur R.
    fprintf(fid, '\n');
  endif
  if out.T
    save_vector(fid, results.T, "T"); % Sauvegarde du vecteur T.
    fprintf(fid, '\n');
  endif
  if out.scattering_matrix
    fprintf(fid, 'Scattering Matrix = \n');
    for i = 1:size(results.S, 1)
    fprintf(fid, '%.4f\t', results.S(i, :));  % Sauvegarde ligne par ligne avec tabulations (\t).
    fprintf(fid, '\n');  % Saut de ligne à la fin de chaque ligne.
    end
    fprintf(fid, '\n'); % Ligne vide entre les vecteurs.

  endif
  if out.conductivity
    fprintf(fid, 'Conductivity = %.4e\n', results.sigma_g); % Sauvegarde de la conductivité.
    fprintf(fid, '\n');  % Saut de ligne à la fin de chaque ligne.

  endif

  fclose(fid); % Fermeture du fichier.

  endif
  return  % Fin de la fonction principale.

  function S = compute_scattering_matrix (gamma1, gamma2, omega, M) % Calcul de la matrice S à l'ordre M de Fourier à partir des matrices gamma1, gamma2, et omega.
  A=eye(4*M+2); % Matrice Identité de taille 4*M+2 car R et T de taille 2*M+1 (de -M à +M en passant par 0).
  B=-A; % Opposé de A ( cf. page 23/24 du diaporama).

  C=gamma1;
  D=gamma2+omega;
  E=gamma2-omega;
  % Conforme à la conforme à la construction de la Matrice S (page 23/24 du diaporama).
  M1=[A, B; C, D];
  M2=[B, A; C, E];
  S=M1\M2; % Algo "division" de deux matrices mieux optimisé que celui de l'inversion.

endfunction

function A = compute_absorption (R, T, M, alpha, beta, d, k1, k2) % Calcul de l'absorption à partir des vecteurs complexes  R et T à l'ordre de troncature de Fourier M. Des composantes de k : alpha, et beta.
% Période d du réseau (m). Des vecteurs d'onde dans le milieu 1 et 2 (k1 et k2).

  res_sum=0;
  gamma_10=sqrt(k1^2 - alpha^2 - beta^2);  % Calcul de gamma_10 : composante z du vecteur d'onde dans le milieu incident
  K=2*pi/d; % Calcul du vecteur d'onde.

  n = (-M:M); % n : indice du canal

  % Indices utilisés pour extraire les composantes x, y, z dans R et T (car sigma dépend de p et pas de n).
  px=n+M+1; % px, py, pz : indices des composantes x, y et z associées au canal n dans le vecteur correspondant
  py=px+1+2*M;
  pz=py+1+2*M;
  alpha_n = transpose(alpha+n*K);
  % Composantes z des vecteurs d'onde diffractés dans les milieux 1 et 2.
  gamma_1n = sqrt(k1^2 - alpha_n.^2 - beta^2);
  gamma_2n = sqrt(k2^2 - alpha_n.^2 - beta^2);
   % Contribution des ondes réflechies et transmises.
  epsilon_rn=real(gamma_1n./gamma_10).*(abs(R(px)).^2 + abs(R(py)).^2 + abs(R(pz)).^2);
  epsilon_tn=real(gamma_2n./gamma_10).*(abs(T(px)).^2 + abs(T(py)).^2 + abs(T(pz)).^2);
  res_sum = sum(epsilon_rn) + sum(epsilon_tn);

  A=1-res_sum; %Calcul de l'absoprtion.
endfunction

function T = toeplitz_matrix (M, a, d, sigma_g) % Fonction qui calcule la matrice de Toeplitz complexe de taille (2M+1) x (2M+1) [de 0 à +-2M).
% M: nombre de strips de graphène . % a : périodicté du réseau. % d période du réseau.
% Fonction sigma_g  : conductivité du graphène.

  % Initialisation :
  c_vect=[]; % Première colonne de la matrice Toeplitz.
  r_vect=[]; % Première ligne de la matrice Toeplitz.
  n = 0:2*M; %Indice qui va de 0 à 2M.


  % Chaque diagonale est constante, définie à partir de c_vect et r_vect
  c_vect = get_sigma_coeff(n, a, d, sigma_g);
  r_vect = get_sigma_coeff(-n, a, d, sigma_g);

  T = toeplitz(c_vect,r_vect);

endfunction

function sigma_g = sigma(w, damp, T, mu_c) % Fonction qui calcule la conductivité de surface du graphène.
% w : la pulsation temporelle (rad/s), damp : l'amortissement (1/s), T (en Kelvin) : la température et le potentiel chimique mu_c (en Joule).

  hbar=1.05E-34; %Constante de Planck réduite (J·s).
  kb=1.38E-23; % constante de Boltzmann (J/K).
  e=1.6E-19; %charge élémentaire (Coulomb)


  sigma_intra = ((2*kb*T*e^2)*i)*(log(2*cosh(mu_c/(2*kb*T))))/(pi*hbar^2*(w + i*damp));
  % Conductivité intrabande dominante à basse fréquence.
  sigma_inter = (e^2/(4*hbar))*(0.5 + atan((hbar*(w + damp*i) -2*mu_c)/(2*kb*T))/pi - i*log(((hbar*(w+i*damp)+2*mu_c)^2)/((hbar*(w+i*damp)-2*mu_c)^2 + (2*kb*T)^2))/(2*pi));
  % Conductivité interbande, dominante à haute fréquence.
  sigma_g = sigma_inter + sigma_intra;
  %Conductivité totale.
endfunction

function [omega, gamma1, gamma2] = omega_gamma_matrix (M, alpha, beta, d, k0, k1, k2, Z0, sigma_g, a) % Fonction qui construit les matrices omega, gamma1 et gamma2.
% Ordre de troncature de Fourier M. Première et deuxième composantes (alpha, beta) du vecteur d'onde. Période d du réseau (m). Différents vecteurs d'ondes (vide : k0, milieu 1 : k1, milieu 2 : k2).

  K=2*pi/d; % Vecteur réciproque liée à la périodicité spatiale.

  Z=zeros(2*M+1); % Création d'une matrice nulle de taille (2M+1)*(2M+1).
  P=k0*Z0*toeplitz_matrix(M, a, d, sigma_g);

  omega=[Z, P; P, Z]; % Assemblage de la matrice omega comme une matrice bloc 2x2.

  A=zeros(2*M+1);
  B=zeros(2*M+1);
  C=zeros(2*M+1);
  D=zeros(2*M+1);
  E=zeros(2*M+1);
  F=zeros(2*M+1);

  for p = 1:2*M+1
    n=p-(M+1); #p matrix index and n true channel index
    alpha_n=alpha+n*K;   %calcule du vecteur alpha_n.
    	%Calcul des gamma des deux milieux.
    gamma_1n=sqrt(k1^2 - alpha_n^2 - beta^2);
    gamma_2n=sqrt(k2^2 - alpha_n^2 - beta^2);

   % Calcul des sous-matrices de gamma1 et gamma2, on exprime chaque sous-bloc comme un produit avec la matrice identité pour traiter les termes diagonaux.
    % Si gamma_1n ou gamma_2n égale à 0, on ajout 0 aux matrices, sinon on divise par zéro
    if gamma_1n == 0
      A(p,p) = 0;
      B(p,p) = 0;
      C(p,p) = 0;
    else
      A(p,p) = alpha_n*beta/gamma_1n;
      B(p, p) = gamma_1n + (beta^2/gamma_1n);
      C(p, p) = gamma_1n + (alpha_n^2/gamma_1n);
    endif

    if gamma_2n == 0
      D(p,p) = 0;
      E(p,p) = 0;
      F(p,p) = 0;
    else
      D(p, p) = alpha_n*beta/gamma_2n;
      E(p, p) = gamma_2n + (beta^2/gamma_2n);
      F(p, p) = gamma_2n + (alpha_n^2/gamma_2n);
    endif
  endfor

  % Assemblage des matrices gamma1 (milieu 1), et gamma2 (milieu2).
  gamma1=[A, B; C, A];
  gamma2=[D, E; F, D];

endfunction

function [I, J, R, T] = IJRT_vectors (Ixy, Jxy, S, M, d, alpha, beta, k1, k2) % Fonction qui permet le calcule des vecteurs complexes I, J, R, T.
% Ixy, Jxy : vecteurs colonnes représentant les composantes  x, y. % S : scattering matrix reliant les vecteurs incidents aux vecteurs réfléchis et transmis.
% M : ordre de troncature pour les harmoniques de Fourier. % d : Période d du réseau (m). % alpha  et beta : composante x  et y du vecteur d’onde incident.
% k1, k2   : modules du vecteur d’onde dans les milieux 1 (incident) et 2 (transmis)

  K=2*pi/d;  % Vecteur réciproque liée à la périodicité spatiale.


  I_J = [Ixy;Jxy]; % Vecteur qui contient les composantes en x et en y de I et J.
  R_T = S*I_J; % Vecteur qui contient les composantes en x et en y de R et T (obtenu comme le produit de la matrice de scattering S par I_J).

  R=R_T(1:4*M+2); % On va de 1 à 4*M+2 car R contiennent 2M+1 valeurs (Rx, Ry).
  T=R_T(4*M+3:8*M+4); % On va de 4*M+3 à 8*M+4 car T contiennent 2M+1 valeurs (Tx, Ty). De plus on reprend l'indice de sommation de R_T incrementé de 1.

  %On degage les composantes en x et y de I et J.
  I = Ixy;
  J = Jxy;

  for px = 1:2*M+1 %px: x vector index
    py = px+1+2*M;  %py: y vector index
    n = px-1-M; %n: channel index

    alpha_n=alpha+n*K; % calcule du vecteur alpha_n.
      %Calcul des gamma des deux milieux.
    gamma_1n=sqrt(k1^2 - alpha_n^2 - beta^2);
    gamma_2n=sqrt(k2^2 - alpha_n^2 - beta^2);

      % Calcul des composantes en z (par la formule page 20/24 du diaporama).
    if gamma_1n == 0
      I(end+1)=0;
      R(end+1)=0;
    else
      I(end+1)=-(alpha_n*I(px) + beta*I(py))/gamma_1n;
      R(end+1)=(alpha_n*R(px) + beta*R(py))/gamma_1n;
    endif

    if gamma_2n == 0
      J(end+1)=0;
      T(end+1)=0;
    else
      J(end+1)=(alpha_n*J(px) + beta*J(py))/gamma_2n;
      T(end+1)=-(alpha_n*T(px) + beta*T(py))/gamma_2n;
    endif

  endfor

endfunction

function sigma_p = get_sigma_coeff (p, a, d, sigma_g) % Calcule les coefficients de la matrice de Toeplitz (ou de convolution) pour une bande de graphène qui est périodique associés aux harmoniques p.
% p : vecteur des indices d’harmoniques de Fourier (entiers). % a : largeur de la bande de graphène. % d : période d du réseau. % sigma_g  : conductivité du graphène.



  sigma_p = []; % Initialisation du vecteur résultat.

  for h = 0:columns(p)-1 % Parcours des indices de Fourier (de 0 à 2M).

    % Harmonique d'ordre 0 : Formule page 12/24.
    if (h==0)
      o = sigma_g * a / d;
	% Coefficient issu de l'intégrale de Fourier : Formule page 12/24.
    else
      K = 2*pi/d;
      o = (exp(-i*K*p(h+1)*a) - 1)*( i*sigma_g/(2*p(h+1)*pi) );
    endif
    sigma_p(end+1) = o; %On remplit le vecteur sigma_p.
  endfor

endfunction

function save_vector(fid, vector, vector_name) % Fonction qui enregistre un vecteur complexe dans un fichier texte ouvert sous forme algebrique.
% fid : identifiant du fichier % vector : vecteur complexe à sauvegarder % vector_name : étiquette du vecteur à écrire dans le fichier.

  fprintf(fid, '\n%s\n', vector_name);  % Écrit le nom du vecteur dans le fichier, suivi d’un saut de ligne.

  data = [real(vector)(:), imag(vector)(:)]';
  fprintf(fid, '%+e%+ei\n', data); % Ecriture du vecteur sous forme algébrique (a+ ib).
endfunction


endfunction % Fin fonction run_simulation ()
