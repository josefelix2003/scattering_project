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

Iz=-(alpha*Ix+beta*Iy)/(gamma1)
Jz=-(alpha*Jx+beta*Jy)/(gamma2)
Tz=-(alpha*Tx+beta*Ty)/(gamma2)
Rz=-(alpha*Rx+beta*Ry)/(gamma1)


% Création des vecteurs I, R et k
I = [Ix; Iy,Iz];
R = [Rx; Ry,Rz];

% Afficher les résultats
disp('Vecteur I:');
disp(I);

disp('Vecteur R:');
disp(R);
