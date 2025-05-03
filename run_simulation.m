
function results = run_simulation ()
  config

  function S = compute_scattering_matrix (gamma1, gamma2, omega, M)
  A=eye(4*M+2);
  B=-A;

  C=gamma1;
  D=gamma2+omega;
  E=gamma2-omega;

  M1=[A, B; C, D];
  M2=[B, A; C, E];

  S=M1\M2;

endfunction

function A = compute_absorption (R, T, M, alpha, beta, d, k1, k2)
  res_sum=0;
  gamma_10=sqrt(k1^2 - alpha^2 - beta^2);
  K=2*pi/d;

  for n = -M:M
    px=n+M+1;
    py=px+1+2*M;
    pz=py+1+2*M;

    alpha_n=alpha+n*K;
    gamma_1n=sqrt(k1^2 - alpha_n^2 - beta^2);
    gamma_2n=sqrt(k2^2 - alpha_n^2 - beta^2);

    epsilon_rn=real(gamma_1n/gamma_10)*(abs(R(px))^2 + abs(R(py))^2 + abs(R(pz))^2);
    epsilon_tn=real(gamma_2n/gamma_10)*(abs(T(px))^2 + abs(T(py))^2 + abs(T(pz))^2);

    res_sum+=epsilon_rn+epsilon_tn;

  endfor

  A=1-res_sum;
endfunction

function T = toeplitz_matrix (M, a, d, sigma_g)
  c_vect=[];
  r_vect=[];
  for n = 0:2*M
    c_vect(end+1) = get_sigma_coeff(n, a, d, sigma_g);
    r_vect(end+1) = get_sigma_coeff(-n, a, d, sigma_g);

    endfor


  T = toeplitz(c_vect,r_vect);

endfunction

function sigma_g = sigma(w, damp, T, mu_c)

  e=1.6E-19;
  hbar=1.05E-34;
  kb=1.38E-23;

  sigma_intra = ((2*kb*T*e^2)*i)*(log(2*cosh(mu_c/(2*kb*T))))/(pi*hbar^2*(w + i*damp));
  sigma_inter = (e^2/(4*hbar))*(0.5 + atan((hbar*(w + damp*i) -2*mu_c)/(2*kb*T))/pi - i*log(((hbar*(w+i*damp)+2*mu_c)^2)/((hbar*(w+i*damp)-2*mu_c)^2 + (2*kb*T)^2))/(2*pi));
  sigma_g = sigma_inter + sigma_intra;

endfunction

function [omega, gamma1, gamma2] = omega_gamma_matrix (M, alpha, beta, d, k0, k1, k2, Z0, sigma_g, a)

  K=2*pi/d;

  Z=zeros(2*M+1);
  P=k0*Z0*toeplitz_matrix(M, a, d, sigma_g);

  omega=[Z, P; P, Z];

  A=zeros(2*M+1);
  B=zeros(2*M+1);
  C=zeros(2*M+1);
  D=zeros(2*M+1);
  E=zeros(2*M+1);
  F=zeros(2*M+1);


  for p = 1:2*M+1
    n=p-(M+1); #p matrix index and n true channel index
    alpha_n=alpha+n*K;
    gamma_1n=sqrt(k1^2 - alpha_n^2 - beta^2);
    gamma_2n=sqrt(k2^2 - alpha_n^2 - beta^2);

    A(p, p) = alpha_n*beta/gamma_1n;
    B(p, p) = gamma_1n + (beta^2/gamma_1n);
    C(p, p) = gamma_1n + (alpha_n^2/gamma_1n);

    D(p, p) = alpha_n*beta/gamma_2n;
    E(p, p) = gamma_2n + (beta^2/gamma_2n);
    F(p, p) = gamma_2n + (alpha_n^2/gamma_2n);



  endfor

  gamma1=[A, B; C, A];
  gamma2=[D, E; F, D];

endfunction

function [I, J, R, T] = IJRT_vectors (Ixy, Jxy, S, M, d, alpha, beta, k1, k2)
  K=2*pi/d;

  I_J = [Ixy;Jxy];
  R_T = S*I_J;

  R=R_T(1:4*M+2);
  T=R_T(4*M+3:8*M+4);

  I = Ixy;
  J = Jxy;

  for px = 1:2*M+1 %px: x vector index
    py = px+1+2*M;  %py: y vector index
    n = px-1-M; %n: channel index

    alpha_n=alpha+n*K;
    gamma_1n=sqrt(k1^2 - alpha_n^2 - beta^2);
    gamma_2n=sqrt(k2^2 - alpha_n^2 - beta^2);


    I(end+1)=-(alpha_n*I(px) + beta*I(py))/gamma_1n;
    J(end+1)=(alpha_n*J(px) + beta*J(py))/gamma_2n;
    R(end+1)=(alpha_n*R(px) + beta*R(py))/gamma_1n;
    T(end+1)=-(alpha_n*T(px) + beta*T(py))/gamma_2n;


  endfor


endfunction

function sigma_p = get_sigma_coeff (p, a, d, sigma_g)
  if (p==0)
    sigma_p = sigma_g * a / d;
  else
    K = 2*pi/d;
    sigma_p = (exp(-i*K*p*a) - 1)*( i*sigma_g/(2*p*pi) );
  endif

endfunction

function results = compute_values (in, out)

  #--------------------Constants-----------------------
  c=3E8;
  hbar=1.05E-34;
  permittivity0=8.8541E-12;
  permeability0=1.2566E-6;
  Z0=sqrt(permeability0/permittivity0);
  damping_factor=0.658E-22/hbar; %à vérifier si c'est une constante
  #----------------------------------------------------
  #-----------------Derived quantities-----------------
  k0=2*pi/in.lambda;
  w=2*pi*c/(in.n*in.lambda); 
  sigma_g=sigma(w, damping_factor, in.T, in.mu_c);
  k1=k0*sqrt(in.permittivityr1*in.permeabilityr1);
  k2=k0*sqrt(in.permittivityr2*in.permeabilityr2);

                  %Incident wave vector
  alpha=k1*sin(in.theta)*cos(in.phi);
  beta=k1*sin(in.theta)*sin(in.phi);
  gamma=k1*cos(in.theta);
  #----------------------------------------------------
  #-------------Computing Scattering Matrix------------
              %Omega, Gamma1, Gamma2 Matrices
  [omega, gamma1, gamma2] = omega_gamma_matrix(in.M, alpha, beta, in.d, k0, k1, k2, Z0, sigma_g, in.a);
                    %Scattering Matrix
  S=compute_scattering_matrix(gamma1, gamma2, omega, in.M);
  #----------------------------------------------------
  #-----------------Compute R, T vectors---------------
  [I, J, R, T]=IJRT_vectors(in.Ixy, in.Jxy, S, in.M, in.d, alpha, beta, k1, k2);
  #----------------------------------------------------
  #--------------------Absorption----------------------
  A=compute_absorption(R, T, in.M, alpha, beta, in.d, k1, k2);
  #----------------------------------------------------
  #-------------Building Results structure-------------
  results=struct();

  if out.absorption
    results.absorption = A;
  endif
  if out.I
    results.I = I;
  endif
  if out.J
    results.J = J;
  endif
  if out.R
    results.R = R;
  endif
  if out.T
    results.T = T;
  endif
  if out.scattering_matrix
    results.scattering_matrix = S;
  endif
  if out.conductivity
    results.conductivity = sigma_g;
  endif

endfunction

  if sweep.enable
    sweep_var=sweep.input_variable;
    output_var=sweep.output_variable;
    min_val=sweep.min;
    max_val=sweep.max;
    step_val=sweep.step;
    xaxis=zeros(0,1);
    yaxis=zeros(0,1);
    counter=0;

    for current_val = min_val:step_val:max_val
      counter+=1;
      in.(sweep_var)=current_val;
      sim_results = compute_values(in, out);
      xaxis(counter,1)= current_val;
      yaxis(counter,1)= sim_results.(output_var);
    endfor
    results=struct();
    results.xaxis=xaxis;
    results.yaxis=yaxis;

  #----------------------------------------------------
  #---------------Building Results structure-----------

   % Obtient le chemin absolu du fichier actuel (script ou fonction)
    current_file_path = mfilename('fullpath');
    [folder_path, ~, ~] = fileparts(current_file_path);

    % Crée le chemin vers le dossier "results" situé dans le même dossier que le script
    result_folder = fullfile(folder_path, 'results');

    % Crée le dossier "results" s’il n’existe pas déjà
    if ~exist(result_folder, 'dir')
        mkdir(result_folder);
    end

    % Spécifie le chemin complet pour enregistrer les résultats dans un fichier texte
    fichier = fullfile(result_folder, 'results.txt');
    fid = fopen(fichier, 'w');

    % Vérifie que le fichier est bien ouvert
    if fid == -1
        error('Impossible d''ouvrir le fichier.');
    end

    % Écrit l'en-tête dans le fichier texte
    fprintf(fid, 'Simulation Results:\n\n');
    fprintf(fid, '%s\t%s\n', sweep_var, output_var);

    % Écrit les résultats ligne par ligne
    for i = 1:length(xaxis)
        fprintf(fid, '%.4e\t%.4e\n', xaxis(i), yaxis(i));
    end

    % Ferme le fichier
    fclose(fid);

    % Si la condition pour afficher le plot est vraie
    if sweep.plot
        % Crée le plot
        plot(xaxis, yaxis);
        xlabel(sweep_var, 'FontSize', 16, 'FontName', 'Times');
        ylabel(output_var, 'FontSize', 16, 'FontName', 'Times');
        grid on;
        set(gca, 'FontSize', 12, 'FontName', 'Times');
        axis tight;
        box on;

        % Génère le titre en fonction des variables
        title_str = sprintf('%s vs. %s', output_var, sweep_var);
        title(title_str, 'FontSize', 16, 'FontWeight', 'bold');

        % Définit le chemin du fichier pour enregistrer la figure
        plot_filename = fullfile(result_folder, 'plot.png');

        % Sauvegarde la figure en PNG
        print(gcf, plot_filename, '-dpng');
    end


    return

  #----------------------------------------------------

  endif

  if ~sweep.enable

  results=compute_values(in, out);

  current_file_path = mfilename('fullpath');
  [folder_path, ~, ~] = fileparts(current_file_path);

    % Crée le chemin vers le dossier "resultats"
  result_folder = fullfile(folder_path, 'results');

    % Crée le dossier s’il n’existe pas
  if ~exist(result_folder, 'dir')
    mkdir(result_folder);
  end

  fichier = fullfile(result_folder, 'results.txt');
  fid = fopen(fichier, 'w');

  if fid == -1
    error('Impossible d''ouvrir le fichier.');
  end

  fprintf(fid, 'Simulation Results:\n\n');

  if out.absorption
    fprintf(fid, 'Absorption = %.4e\n', results.absorption);
  endif
  if out.I
    fprintf(fid, 'I = \n');
    for i = 1:size(results.I, 1)
    fprintf(fid, '%.4f\t', results.I(i, :));  % "\t" pour tabulation entre les valeurs
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne
    end
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne

  endif
  if out.J
    fprintf(fid, 'J = \n');
    for i = 1:size(results.J, 1)
    fprintf(fid, '%.4f\t', results.J(i, :));  % "\t" pour tabulation entre les valeurs
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne
    end
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne

  endif
  if out.R
    fprintf(fid, 'R = \n');
    for i = 1:size(results.R, 1)
    fprintf(fid, '%.4f\t', results.R(i, :));  % "\t" pour tabulation entre les valeurs
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne
    end
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne

  endif
  if out.T
    fprintf(fid, 'T = \n');
    for i = 1:size(results.T, 1)
    fprintf(fid, '%.4f\t', results.T(i, :));  % "\t" pour tabulation entre les valeurs
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne
    end
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne

  endif
  if out.scattering_matrix
    fprintf(fid, 'Scattering Matrix = \n');
    for i = 1:size(results.S, 1)
    fprintf(fid, '%.4f\t', results.S(i, :));  % "\t" pour tabulation entre les valeurs
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne
    end
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne

  endif
  if out.conductivity
    fprintf(fid, 'Conductivity = %.4e\n', results.sigma_g);
    fprintf(fid, '\n');  % Nouvelle ligne à la fin de chaque ligne

  endif

  fclose(fid);

  endif
  return
endfunction



