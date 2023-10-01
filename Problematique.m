clc;
close;
clear all;

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

%% Matrice ABCD FTBO
A = [[0 0 0 0]' [1 -(Bm*N^2*Bl)/(N^2*Jl+Jm) Kb/(N*La) 0]' [0 (N*Ki)/(N^2*Jl+Jm) -Ra/La 0]' [0 0 1/La -1/tau]'];
B = [ 0 0 0 K/tau]';
C = [1 0 0 0];
D = [0];
[num,denum] = ss2tf(A,B,C,D)
FTBO = tf(num,denum)

%% Matrice ABCD FTBF
A_FBTF = A;
A_FBTF(4,:) = [(-K*Kp)/tau 0 0 -1/tau]';

B_FBTF = [0 0 0 (K*Kp)/tau]';
[num_FBTF,denum_FBTF] = ss2tf(A_FBTF, B_FBTF, C, D)
FBTF = tf(num_FBTF,denum_FBTF)
 
