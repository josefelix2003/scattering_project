in=struct();
out=struct();
sweep=struct();
#-------------------Sweep config----------------------
sweep.enable = true; % Active/désactive le balayage
sweep.input_variable = 'lambda'; % Nom de la variable d’entrée à balayer
sweep.min = 300E-9;
sweep.max = 400E-9;
sweep.step = 50E-9;
sweep.output_variable = 'absorption'; % Nom de la grandeur à extraire
sweep.plot=true;

#-----------------------------------------------------
#-----------------Input parameters--------------------
                %Grating parameters
in.M=4;
in.a=1E-9;
in.d=3E-9;
             %Relative Material parameters
in.n=1; %à vérifier
in.permittivityr1=1;
in.permeabilityr1=1;
in.permittivityr2=1;
in.permeabilityr2=1;
in.T=300;
in.mu_c=5E-20; %en Joules
               %Incident Wave parameters
in.delta=pi/3;
in.phi=pi/4;
in.theta=pi/6;
in.lambda=400E-9;
              %I, J incident vectors
in.I0_x=cos(in.delta)*cos(in.phi)*cos(in.theta)-sin(in.delta)*sin(in.phi);
in.I0_y=cos(in.delta)*sin(in.phi)*cos(in.theta)+sin(in.delta)*cos(in.phi);
in.Ixy=zeros(4*in.M+2,1);
in.Jxy=zeros(4*in.M+2,1);
in.Ixy(in.M+1)=in.I0_x;
in.Ixy(3*in.M+2)=in.I0_y;
#-----------------------------------------------------
#-----------------Output parameters-------------------
out.absorption=true;
out.I=true;
out.J=false;
out.R=false;
out.T=false;
out.scattering_matrix=false;
out.conductivity=false;


