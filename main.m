#Input parameters
M=3;
a=0.0001;
d=0.1;
delta=pi/2;
phi=pi/2;
theta=pi/2;
lambda=
k0= 2*pi/lambda %vérifier si c'est une inconnue
permittivity1=
permeability1=
permittivity2=
permeability2=

permittivity0=8.8541E-12
permeability0=1.2566E-6
Z0=sqrt(permeability0/permittivity0)

sigma_g=1.5; %fonction pour calculer la valeur de sigma_g

#Incident wave vector
alpha=sin(theta)*cos(phi)
beta=sin(theta)*sin(phi)
gamma=cos(theta)

k1=k0*sqrt(permittivity1*permeability1)
k2=k0*sqrt(permittivity2*permeability2)



#I0
I0=[cos(delta)*cos(theta)*cos(phi) - sin(delta)*sin(phi); cos(delta)*cos(theta)*sin(phi) + sin(delta)*sin(phi); -cos(delta)*sin(theta)];

#I and J vectors
J=zeros(6*M+3, 1);
I=zeros(6*M+3, 1);
I(M+1)=I0(1);
I(3*M+2)=I0(2);
I(5*M+3)=I0(3);




#Toeplitz matrix

T = toeplitz(c_vect,r_vect)
