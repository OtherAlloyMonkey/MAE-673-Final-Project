clc; clear all; close all;

m = 1;
ms = 0.1;
k = 1;
ks = 50;

Mmat = [m+ms ms; ms ms];
kmat = [k 0; 0 ks];
Dmat = [k;0];

A = [zeros(2) eye(2); -inv(Mmat)*kmat zeros(2)];
B = [0;0;inv(Mmat)*Dmat];
C = [eye(2) zeros(2)];
D = [0;0];
eA = eig(A);
sys = ss(A,B,C,D);

% Designing a Time Delay Filter
om1 = abs(eA(end));
% because system is undamped 
T = pi/om1;
A0 = 0.5;
A1 = 0.5;

tvec = linspace(0,5*T,5001);
uvec = ones(size(tvec)); % step input
uvecF = 0.5*ones(size(tvec)) + 0.5*heaviside(tvec-T); % filtered version
figure(1)
plot(tvec,uvec,tvec,uvecF)

[yy, tt, xx] = lsim(sys,uvec,tvec);
[yyF, ttF, xxF] = lsim(sys,uvecF,tvec);

figure(2)
subplot(211),plot(tt,yy(:,1),ttF,yyF(:,1)); hold on
subplot(212),plot(tt,yy(:,2),ttF,yyF(:,2)); hold on

hold on
J = 5;
T1 = 1/(2*J);
T2 = pi/(2*om1) + 1/(4*J);
% creating u-derivative vec 
uDvec = J*ones(size(tvec)) - J*heaviside(tvec-T1) + J*heaviside(tvec-(2*T2-T1)) - J*heaviside(tvec-2*T2);

figure(3)
plot(tvec,uDvec)

sysf = tf([1],[1 0]); %TF of integrator 
uTDFJ = lsim(sysf,uDvec,tvec);  % jerk limimted shaped input
plot(tvec,uTDFJ)

[yyF2, ttF2, xxF2] = lsim(sys,uTDFJ,tvec);
figure(2)
subplot(211),plot(ttF2,yyF2(:,1),'k--')
subplot(212),plot(ttF2,yyF2(:,2),'k--')
