%% Swenson MAE 673 HW1.2 
clear; close all; clc; format compact; format long;

%% Design the state feedback control for the system

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

[num, den] = ss2tf(A,B,C,D);


%% Design TDF for each mode of the OL system

% Need the natrural frequencies of the OL system
sysOL = ss(A,B,C,D); 
wn = imag(pole(sysOL));

T1 = pi/wn(1);
T2 = pi/wn(3);

t = 0:.001:(T1+T2+5);


sysstep = step(sysOL,t);
utdf1 = .5*(1 + heaviside(t-T1));
utdf2 = .5*(1 + heaviside(t-T2));
utdfcasc = den(5)*.25*(1 + heaviside(t-T1) + heaviside(t-T2) + heaviside(t-(T1+T2)) );

xtdf1 = lsim(sysOL,utdf1,t);
xtdf2 = lsim(sysOL,utdf2,t);
[xtdfcasc, tout, state] = lsim(sysOL,utdfcasc,t);


%% Concurrent Design of the TDF to minimize T2

% Natural Frequencies of the closed loop TF

polesCL = pole(SSsys);

% Inequality constraints
Aineq = [];
Bineq = [];

% Equality Constraints
Aeq = [1 1 1 0 0];
Beq = 1;

% Lower and upper bounds
LB = []';
UB = []';

% Function to be minimized (T2) and nonlinear constraint as function of
% just x.
optfun = @(x) x(5);
nonlcon = @(x) nonlincon(x,polesCL); 

% Initial Conditions
x0 = [1/3 1/3 1/3 T1 T2];

% Optimize A_i and T_i
xout = fmincon(optfun,x0,Aineq,Bineq,Aeq,Beq,LB,UB,nonlcon);
A0 = xout(1); A1 = xout(2); A2 = xout(3); T1c = xout(4);   T2c = xout(5);

% Simulate response to concurrent TDF
ucon = den(5)*(xout(1) + xout(2)*heaviside(t-xout(4)) + xout(3)*heaviside(t-xout(5)));
[xconc, tout, statecon] = lsim(SSsys,ucon,t);




figure();
plot(t,xtdfcasc)
hold on
plot(t,xconc)
legend('Cascaded TDF','Concurrent TDF','location','southeast')
title('Displacemnt of Mass subject to TDF')
xlabel('Time in (s)')
ylabel('Position (m)')

%% Create added robustness by cascading the concurrent TDF
A0 = xout(1); A1 = xout(2); A2 = xout(3);
T1c = xout(4);   T2c = xout(5);

t = 0:.1:(2*T2c+10);


ucon = den(5)*(xout(1) + xout(2)*heaviside(t-xout(4)) + xout(3)*heaviside(t-xout(5)));

uconcasc2 = den(5)*((A0^2) + (A1^2)*heaviside(t-2*T1c) + (A2^2)*heaviside(t-2*T2c) +...
            2*A0*A1*heaviside(t-T1c) + 2*A0*A2*heaviside(t-T2c) + 2*A1*A2*heaviside(t-T1c-T2c));

uconcasc3 = den(5)*(A0^3 + heaviside(t-3*T1c)*A1^3 + heaviside(t-3*T2c)*A2^3 + 3*(A0^2)*A1*heaviside(t-T1c) + ...
                    3*(A0^2)*A2*heaviside(t-T2c) + 3*A0*(A1^2)*heaviside(t-2*T1c) + 3*A0*(A2^2)*heaviside(t-2*T2c)...
                  + 3*(A1^2)*A2*heaviside(t-(2*T1c+T2c)) + 3*A1*(A2^2)*heaviside(t-(T1c+2*T2c)) + 6*A0*A1*A2*heaviside(t-T1c-T2c));


k = linspace(.8,1.2,101);
Cost = zeros(length(k),1);
Costc2 = zeros(length(k),1);
Costc3 = zeros(length(k),1);


for ii = 1:length(k)

    % Define a new state space A mat, design new poles for it.
    Atemp = [0 1 0 0; (-2*k(ii))/m 0 r/m 0; 0 0 0 1; 1/J 0 (-k(ii)*r^2)/J 0];
   
    
    system = ss(Atemp,B,C,D);
     
    [ytempcon, ttempcon, statetempcon] = lsim(system,ucon,t);
    [ytempcasc, ttempcasc, statetempcasc] = lsim(system,uconcasc2,t);
    [ytempcasc2, ttempcasc2, statetempcasc3] = lsim(system,uconcasc3,t);

    Cost(ii) = .5*m*(statetempcon(end,2))^2 + .5*J*statetempcon(end,4)^2 + .5*k(ii)*(statetempcon(end,1)-1)^2 + ...
              (.5*k(ii)*(statetempcon(end,1)-1)^2 - .5*k(ii)*(r*statetempcon(end,3) - statetempcon(end,1)- 1)^2 )^2;

    Costc2(ii) = .5*m*(statetempcasc(end,2))^2 + .5*J*statetempcasc(end,4)^2 + .5*k(ii)*(statetempcasc(end,1)-1)^2 + ...
                (.5*k(ii)*(statetempcasc(end,1)-1)^2 - .5*k(ii)*(r*statetempcasc(end,3) - statetempcasc(end,1) - 1)^2 )^2;

    Costc3(ii) = .5*m*(statetempcasc3(end,2))^2 + .5*J*statetempcasc3(end,4)^2 + .5*k(ii)*(statetempcasc3(end,1)-1)^2 + ...
                (.5*k(ii)*(statetempcasc3(end,1)-1)^2 - .5*k(ii)*(r*statetempcasc3(end,3) - statetempcasc3(end,1) - 1)^2 )^2;

end

Cost = Cost.^(1);
Costc2 = Costc2.^(1);
Costc3 = Costc3.^(1);

figure();
semilogy(k,Cost,'o',k,Costc2,'.',k,Costc3,'.')
legend('Concurrent Filter','2 Concurrent Filter Cascade ','3 Concurrent Filter Cascade','location','northeast')
xlabel('Spring stiffness')
ylabel('Square Root Residual Energy Magnitude','interpreter','latex')
title('Filter Sensitivity to Spring Stiffness Uncertainty')
%% Functions

% Function of nonlinear constraints
function [c,ceq] = nonlincon(x,polesCL)

    sigs = -1*real(polesCL);
    wd = imag(polesCL);


    c = [];

    % Nonlinear constraints for CL mode 1
    ceq(1,1) = x(1) + x(2)*exp(sigs(1)*x(4))*cos(wd(1)*x(4)) + x(3)*exp(sigs(1)*x(5))*cos(wd(1)*x(5));
    ceq(2,1) = x(2)*exp(sigs(1)*x(4))*sin(wd(1)*x(4)) + x(3)*exp(sigs(1)*x(5))*sin(wd(1)*x(5));

    % Nonlinear constraints for CL mode 1
    ceq(3,1) = x(1) + x(2)*exp(sigs(3)*x(4))*cos(wd(3)*x(4)) + x(3)*exp(sigs(3)*x(5))*cos(wd(3)*x(5));
    ceq(4,1) = x(2)*exp(sigs(3)*x(4))*sin(wd(3)*x(4)) + x(3)*exp(sigs(3)*x(5))*sin(wd(3)*x(5));

end
