% Demander à l'utilisateur d'entrer les valeurs
Ix = input('Veuillez entrer Ix: ');
Ix = double(Ix);  % Conversion en type flottant

Iy = input('Veuillez entrer Iy: ');
Iy = double(Iy);  % Conversion en type flottant

Rx = input('Veuillez entrer Rx: ');
Rx = double(Rx);  % Conversion en type flottant

Ry = input('Veuillez entrer Ry: ');
Ry = double(Ry);  % Conversion en type flottant

beta = input('Veuillez entrer beta: ');
beta = double(beta);  % Conversion en type flottant

% Création des vecteurs I, R et k
I = [Ix; Iy];
R = [Rx; Ry];
k = [alpha_n; beta; gamma];

% Afficher les résultats
disp('Vecteur I:');
disp(I);

disp('Vecteur R:');
disp(R);

disp('Vecteur k:');
disp(k);
