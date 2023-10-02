clc;
clear;
close all;


%% constant
Kp  = 0.318; %V/rad
K   = 100;
tau = 0.01; %s
Ki  = 0.5; % N*m/A
Kb  = 0.5; %V/rad/s
Ra  = 8; %ohm
La  = 0.008; % H
Jm  = 0.02; %N*m*s^2/rad
Bm  = 0.01; %N*m*s^2/rad
N   = 0.1;
Jl  = 1; %N*m*s^2/rad
Bl  = 1; %N*m*s^2/rad

%% c Matrice ABCD FTBO
A = [[0 0 0 0]' [1 (-N*Bl-Bm/N)/(Jm/N+Jl*N) -Kb/(N*La) 0]' [0 Ki/(Jm/N+Jl*N) -Ra/La 0]' [0 0 1/La -1/tau]'];
B = [ 0 0 0 K/tau]';
C = [1 0 0 0];
D = [0];

%% e fonction de transfert en boucle ouverte 
[num,denum] = ss2tf(A,B,C,D);
FTBO = tf(num,denum)
figure;
pzmap(A,B,C,D);
figure;
impulse(ss(A,B,C,D))

%% Matrice ABCD FTBF
A_FTBF = A;
A_FTBF(4,:) = [(-K*Kp)/tau 0 0 -1/tau]';
B_FTBF = [0 0 0 (K*Kp)/tau]'; 

%% e fonction de transfert en boucle fermee
[num_FTBF,denum_FTBF] = ss2tf(A_FTBF, B_FTBF, C, D);
FTBF = tf(num_FTBF,denum_FTBF)

%% f1 Reduction Physique


%% f2 Reduction Numerique

%reduction a une 2e ordre
[R,P,K] = residue(num,denum) %utiliser FTBO le gain DC est infini
ratio = abs((R)./real(P))
[numr, denumr] = residue(R(3:4),P(3:4),K) % prenre lui infini
%numr = numr*(dcgain(FTBO)/dcgain(num_FTBO,denum_FTBO)) %negliger etape de dcgain
figure;
FTBOr = tf(num,denum)
step(FTBOr,[0:1/1000:10]),hold;
step(FTBOr,[0:1/1000:10]);


%% g reponse FTBF a un step
figure;
step(num_FTBF,denum_FTBF)

%% h E1
load('donnees_moteur_2016.mat');
acceleration = diff(vitesse)./diff(t);
d3 = diff(acceleration)./diff(t(1:4000));
mX = [d3, acceleration(1:3999), vitesse(1:3999) ];
out = pinv(mX)*tension(1:3999)
figure;
plot(t(1:4000),vitesse(1:4000)), hold on;
plot(t(1:4000), acceleration);
%% h E2
V = 8; %V
Tm = 0.52; %Nm
Ia = 1.09; %A
