%% Swenson MAE 673 HW 1
clear; close all; clc; format long;

%% Simulate Open Loop system
J = 1; 
m = 1;
r = 1;
k = 1;

a = J*m/(k*r);
b = 2*J/r + m*r;
c = k*r;



% Create OL Plant TF
Numerator = [1];
Denominator = [a 0 b 0 c];
sys = tf(Numerator,Denominator);
[Aol,Bol,Col,Dol] = tf2ss(Numerator,Denominator);
% SSsys = ss(A,B,C,D);
% stepsim = step(sys);




% Create the TDF Controller TF
g = (3-sqrt(5))/2;  % Slower mode wn ~ .6 rad/s
h = (3+sqrt(5))/2;  % Faster mode wn ~ 1.6 rad/s

T1 = pi/sqrt(g);
T2 = pi/sqrt(h); % Note that I called this T2 but it is the shorter delay

tend = max([T1 T2]) + 10;
t = [0:.05:tend];

u1 = 1*(ones(length(t),1));
u2 = .5*(1 + heaviside(t-T1) );
u3 = .5*(1 + heaviside(t-T2));
u4 = .25*(1+ heaviside(t-T1) + heaviside(t-T2) + heaviside(t-(T1+T2)) );


% State space rep
A = [0 1 0 0; -(2*k/m) 0 r/m 0; 0 0 0 1; 1/J 0 -(k*r^2)/J 0];
B = [0 0 0 1]';
C = [1 0 0 0];
D = 0;
controllability = rank([B A*B (A^2)*B (A^3)*B]);
SSreal = ss(A,B,C,D);

[yr1, tOutr1, xstateCascr1] = lsim(SSreal,u1,t);
[yr2, tOutr2, xstateCascr2] = lsim(SSreal,u2,t);
[yr3, tOutr3, xstateCascr3] = lsim(SSreal,u3,t);
[yr4, tOutr4, xstateCascr4] = lsim(SSreal,u4,t);

figure();
plot(t,xstateCascr1(:,1));
hold on
plot(t,xstateCascr2(:,1));
plot(t,xstateCascr3(:,1));
plot(t,xstateCascr4(:,1));
legend('Step','Slow Mode Filtered','Fast Mode Filtered','Both Filtered','location','west')


%% Concurrent TDF design

% Inequality constraints
Aineq = [0 0 0 -1 0; 0 0 0 1 -1];
Bineq = [0 0]';

% Equality Constraints
Aeq = [1 1 1 0 0];
Beq = 1;

% Lower and upper bounds
LB = [0 0 0 0 0]';
UB = [1 1 1 inf inf]';

% Function to be minimized (T2) and nonlinear constraint as function of
% just x.
optfun = @(x) x(5);
nonlcon = @(x) nonlincon(x,h,g); 

% Initial Conditions
x0 = [1/3 1/3 1/3 T1 T2];

xout = fmincon(optfun,x0,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon);

% Simulate the opotimized concurrent TDF
ucon = [xout(1) + xout(2)*heaviside(t-xout(4)) + xout(3)*heaviside(t-xout(5))];
[ycon, tcon, xstatecon] = lsim(SSreal,ucon,t);

figure(gcf);
plot(t,xstatecon(:,1),'*');




%% Functions

% Function of nonlinear constraints
function [c,ceq] = nonlincon(x,h,g)

    A0c = x(1);
    A1c = x(2);
    A2c = x(3);
    T1c = x(4);
    T2c = x(5);

    sh = sqrt(h);
    sg = sqrt(g);

    c = [];

%     ceq(1,1) = x(1) + x(2)*cos(sg*x(4)) + x(3)*cos(sh*x(5));
%     ceq(2,1) = x(2)*sin(sg*x(4)) + x(3)*sin(sh*x(5));
%     
%     ceq(1,1) = x(1) + x(2)*( -cos(sh*x(4)) + cos(sg*x(4)) ) + x(3)*( -cos(sh*x(5)) + cos(sg*x(5)) );
%     ceq(2,1) = x(2)*( -sin(sh*x(4)) + sin(sg*x(4)) ) + x(3)*( -sin(sh*x(5)) + sin(sg*x(5)) );

    ceq(1,1) = x(1) + x(2)*( cos(sh*x(4)) ) + x(3)*( cos(sh*x(5)) );
    ceq(2,1) = x(2)*( sin(sh*x(4)) ) + x(3)*( sin(sh*x(5)) );
    ceq(3,1) = x(1) + x(2)*( cos(sg*x(4)) ) + x(3)*( cos(sg*x(5)) );
    ceq(4,1) = x(2)*( sin(sg*x(4)) ) + x(3)*( sin(sg*x(5)) );


end



