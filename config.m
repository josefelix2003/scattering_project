in = struct();       % Structure contenant les paramètres d’entrée.
out = struct();      % Structure définissant les sorties à enregistrer.
sweep = struct();    % Structure pour configurer le balayage paramétrique.

#-------------------Sweep config----------------------
sweep.enable = true;                      % Active (true) ou désactive (false) le balayage paramétrique.
sweep.input_variable = 'lambda';         % Nom de la variable d’entrée à faire varier parmi les éléments de la structure input.
sweep.min = 300E-9;                      % Valeur minimale à faire varier.
sweep.max = 400E-9;                      % Valeur maximale à faire varier.
sweep.step = 50E-9;                      % Pas d’incrément.
sweep.output_variable = 'absorption';   % Nom de l'output.
sweep.plot = true;                      % Si true : affiche un graphique du résultat du balayage.

#-----------------------------------------------------
#-----------------Input parameters--------------------


in.M = 4;                   % Ordre de troncature (nombre de strips de graphène).
in.a = 1E-9;                % Largeur du strip de graphène (mètres).
in.d = 3E-9;                % Période du réseau (mètres).

             % Paramètres relatifs aux milieux :
in.permittivityr1 = 1;      % Permittivité relative du milieu 1 (sans unité).
in.permeabilityr1 = 1;      % Perméabilité relative du milieu 1 (sans unité).
in.permittivityr2 = 1;      % Permittivité relative du milieu 2 (sans unité).
in.permeabilityr2 = 1;      % Perméabilité relative du milieu 2 (sans unité).
in.T = 300;                 % Température (Kelvin).
in.mu_c = 5E-20;            % Potentiel chimique (Joules).

               % Paramètres de l’onde incidente
in.delta = pi/3;            % Angle delta de polarisation (rad).
in.phi = pi/4;              % Angle azimutal (rad).
in.theta = pi/6;            % Colatitude (rad).
in.lambda = 400E-9;         % Longueur d’onde incidente (mètres).

              % Vecteurs I et J (champs incidents) :
in.I0_x = cos(in.delta)*cos(in.phi)*cos(in.theta) - sin(in.delta)*sin(in.phi);  % Composante x de I0 (onde incidente dans le canal central).
in.I0_y = cos(in.delta)*sin(in.phi)*cos(in.theta) + sin(in.delta)*cos(in.phi);  % Composante y de I0 (onde incidente dans le canal central).

      % Initialisation des vecteurs (matrice colonne) Ixy et Jxy.
in.Ixy = zeros(4*in.M+2, 1)     % Vecteur qui contient les composantes x de l'onde I en premier et les composantes y en deuxième (par défaut, zéros partout sauf pour le canal central).
in.Jxy = zeros(4*in.M+2, 1);    % Vecteur qui contient les composantes x de l'onde J en premier et les composantes y en deuxième (par défaut zéros partout).

in.Ixy(in.M+1) = in.I0_x;       % Attribution des composantes de I0 dans le vecteur I.
in.Ixy(3*in.M+2) = in.I0_y;

#-----------------------------------------------------
#-----------------Output parameters-------------------

%Choix des paramètres de sortie à activer :
out.absorption = true; % Absorption (sans unité).
out.I = true;          % Vecteur I.
out.J = false;         % Vecteur J.
out.R = false;         % Vecteur R.
out.T = false;         % Vecteur T.
out.scattering_matrix = false;  % Matrice S.
out.conductivity = false;       % Conductivité (S/m).
