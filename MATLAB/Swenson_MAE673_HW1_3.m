%% MAE 673 HW 1 Rev 3
clear; clc; close all; format compact; format long;

%% Start by creating the models for the OL mass and OL wheel plants

J = 1;  % inertia
m = 1;  % mass
k = 1;  % Spring 
r = 1;  % radius

ThetabyUnum = [m 0 2*k];
ThetabyUden = [J*m 0 (2*J*k+m*k) 0 2*(k*r)^2-k*r^2];

% Take this plant to the PID Tuner and find Kp,Ki,Kd values
ThetaU_OL = tf(ThetabyUnum,ThetabyUden);
C = load('PID_C_HW1.mat');

Kp = C.C.Kp;
Ki = C.C.Ki;
Kd = C.C.Kd;

% Create Closed loop plant
TU_CL_num = [Kd Kp (2*Kd+Ki) 2*Kp 2];
TU_CL_den = [1 Kd (3+Kp) (2*Kd+Ki) (2*Kp + 1) 2*Ki];

TU_CL = tf(TU_CL_num,TU_CL_den);
XT_CL = tf([1],[1 0 2])

G_CL = TU_CL*XT_CL

Fg = .1028;


%% Design TDF for each mode of the OL system

% Physical Parameters of the system
J = 1;      % Moment of inertia of the disc
m = 1;      % Mass of block
r = 1;      % Radius of the disc
k = 1;      % Stiffness of each spring

% State space rep for the physical OL system
A = [0 1 0 0; (-2*k)/m 0 r/m 0; 0 0 0 1; 1/J 0 (-k*r^2)/J 0];
B = [0 0 0 1]';
C = [1 0 0 0];
D = 0;
SSsys = ss(A,B,C,D);

% Need the natrural frequencies of the OL system
sysOL = ss(A,B,C,D); 
wn = imag(pole(sysOL));

T1 = pi/wn(1);
T2 = pi/wn(3);

t = 0:.0001:(T1+T2+16);


sysstep = step(sysOL,t);
utdf1 = .5*(1 + heaviside(t-T1));
utdf2 = .5*(1 + heaviside(t-T2));
utdfcasc = .25*(1 + heaviside(t-T1) + heaviside(t-T2) + heaviside(t-(T1+T2)) );

xtdf1 = lsim(sysOL,utdf1,t);
xtdf2 = lsim(sysOL,utdf2,t);
[xtdfcasc, tout, state] = lsim(sysOL,utdfcasc,t);

[x,t,s] = lsim(sysCL,[1.0301*ones(length(t),1)],t);
figure();
plot(x)