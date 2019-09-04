%par for for 6x intervals and symbolic to global variable
Mz0 = 100; % Mz0 is in function also
%tspan = 0:.01:1; %can be [0 1]
t_0=0.0; t_1=5; 
dt=10^-6;
%dt=0.00001;
tspan = t_0:dt:t_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pB=0.215343;  %see also in function
pF=0.784436;   % n = 0.5
pE=0.000221; 
concentration=0.5;

%  pB=0.107683;  %see also in function
%  pF=0.892206;   % n = 0.25
%  pE=0.000110;
%  concentration=0.25;

% pB=0.04308;  %see also in function
% pF=0.95688;   % n = 0.1
% pE=4.4E-5;
% concentration=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 MxB = 0 ;
 MyB = Mz0*pB; 
 MzB = 0;
 MxF = 0;
 MyF = Mz0*pF;
 MzF = 0;
 MxE = 0;
 MyE = Mz0*pE ;
 MzE = 0;

 U0(1) = MxB;    
 U0(2) = MyB;    
 U0(3) = MzB;    
 U0(4) = MxF;    
 U0(5) = MyF;    
 U0(6) = MzF;    
 U0(7) = MxE;    
 U0(8) = MyE;    
 U0(9) = MzE; 
 
 opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
 [t_out,OUT] = ode45(@SSSolver,tspan,U0, opts);

%%%%%%%%%%%%%%  ????????? ????????

% t0 = 0:.01:0.5; 
% y0 = interp1(t_out,OUT,t0);
% plot(t0,y0(:,1))

%%%%%%%%%%%%%%%

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman');


figure(666) 
suptitle(strcat(['Concentration n =',num2str(concentration),' step = ',num2str(dt)]))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.07, 1, 0.93]);
subplot(2,2,1);
plot(t_out, OUT(:,1), 'r'); % can be |1:size(OUT,1)| as Ox
hold on
plot(t_out, OUT(:,2), 'b');
plot(t_out, OUT(:,3), 'g');
xlabel('time (s)')
ylabel('Magnetization')
legend('M_x B','M_y B','M_z B')

%figure(667)
subplot(2,2,2);
plot(t_out, OUT(:,4), 'r'); % can be |1:size(OUT,1)| as Ox
hold on
plot(t_out, OUT(:,5), 'b');
plot(t_out, OUT(:,6), 'g');
xlabel('time (s)')
ylabel('Magnetization')
legend('M_x F','M_y F','M_z F')

%figure(668) 
subplot(2,2,3);
plot(t_out, OUT(:,7), 'r'); % can be |1:size(OUT,1)| as Ox
hold on
plot(t_out, OUT(:,8), 'b');
plot(t_out, OUT(:,9), 'g');
xlabel('time (s)')
ylabel('Magnetization')
legend('M_x E','M_y E','M_z E')

%figure(669)
subplot(2,2,4);
plot(t_out, OUT(:,1)+OUT(:,4)+OUT(:,7), 'r'); % can be |1:size(OUT,1)| as Ox
hold on
plot(t_out, OUT(:,2)+OUT(:,5)+OUT(:,8), 'b');
plot(t_out, OUT(:,3)+OUT(:,6)+OUT(:,9), 'g');
xlabel('time (s)')
ylabel('Magnetization')
legend('M_x ALL','M_y ALL','M_z ALL')

function [ output_args ] = SSSolver(tt,inp_var)

MxB=inp_var(1);
MyB=inp_var(2);
MzB=inp_var(3);
MxF=inp_var(4);
MyF=inp_var(5);
MzF=inp_var(6);
MxE=inp_var(7);
MyE=inp_var(8);
MzE=inp_var(9);

w1 = 0; % perpendiculiar RF field
%w1 = 400*10^6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pB=0.215343;  %see also in function
pF=0.784436;   % n = 0.5
pE=0.000221;   

% pB=0.107683;  %see also in function
% pF=0.892206;   % n = 0.25
% pE=0.000110;

% pB=0.04308;  %see also in function
% pF=0.95688;   % n = 0.1
% pE=4.4E-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%kBF = 1;
%kFB = 1;
%kBF = 1000; 
%kFB = 1000;

%kBF = 0;
kBF = 10^6;
kFB = kBF*(pB+pE)/pF;

%kEB = 1;
%kBE = 1;

%kEB = 0;
kEB = 17.582084;
kBE = (pE/(pB+pF))*kEB; %

%kEB = 10^-1;
%kBE = 10^-1;

deltaB = 4.7; % 4.7 ppm on 500MHz 2335Hz
deltaF = 4.7;
deltaE = 8.3; %8.3 \ 4150 Hz

R1B = 10; %0.1 s
R2B = 10;

R1F = 0.33; %3s
R2F = 0.33;

R1E = 10;
R2E = 10;

Mz0 = 100;

%syms R2B R1B R2F R1F R2E R1E kEB kBE kFB kBF deltaB deltaF deltaE w1 MxB MyB MzB MxF MyF MzF MxE MyE MzE Mz0

 BM = [-R2B-kBF-kBE -deltaB 0 kFB 0 0 kEB 0 0; deltaB -R2B-kBF-kBE -w1 0 kFB 0 0 kEB 0; 
    0 w1 -R1B-kBF-kBE 0 0 kFB 0 0 kEB; kBF 0 0 -R2F-kFB -deltaF 0 0 0 0; 0 kBF 0 deltaF -R2F-kFB -w1 0 0 0; 0 0 kBF 0 w1 -R1F-kFB 0 0 0;
    kBE 0 0 0 0 0 -R2E-kEB -deltaE 0; 0 kBE 0 0 0 0 deltaE -R2E-kEB -w1; 0 0 kBE 0 0 0 0 w1 -R1E-kEB;];

%M= [MxB;MyB;MzB-Mz0;MxF;MyF;MzF-Mz0;MxE;MyE;MzE-Mz0];  % Case of Mz0 for all equal

M= [MxB;MyB;MzB-Mz0*pB;MxF;MyF;MzF-Mz0*pF;MxE;MyE;MzE-Mz0*pE]; % Case of Mz0 is split equally


%%%%%%%%%%% GPU MULTIPLY
% BMM=gpuArray(BM);
% MM=gpuArray(M);
% x=BMM*MM;
% X=gather(x);
%%%%%%%%%%%%%%%%%%%%

X=BM*M;


%lambda=1/(trace((BM)^-1))

 output_args(1) = X(1);
 output_args(2) = X(2);
 output_args(3) = X(3);
 output_args(4) = X(4);
 output_args(5) = X(5);  
 output_args(6) = X(6);
 output_args(7) = X(7);
 output_args(8) = X(8);
 output_args(9) = X(9);

 output_args=output_args';
end 