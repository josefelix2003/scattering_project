% Demander les valeurs
Ix = input('Veuillez entrer Ix: ');
Ix = double(Ix);  % Conversion en type flottant

Iy = input('Veuillez entrer Iy: ');
Iy = double(Iy); 

Rx = input('Veuillez entrer Rx: ');
Rx = double(Rx); 

Ry = input('Veuillez entrer Ry: ');
Ry = double(Ry); 

beta = input('Veuillez entrer beta: ');
beta = double(beta);  

gamman = zeros(1, 2*M + 1);    
n_val = -M:M;

for i = 1:(2*M + 1)
    n = n_val(i);
    alpha_n = alpha + n * 2 * pi / d;
    gamma_n = sqrt(k_p^2 - alpha_n^2 - beta^2);
    gamman(i) = gamma_n;
end


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
