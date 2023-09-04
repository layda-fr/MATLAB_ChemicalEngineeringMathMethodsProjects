%   Get started
clear all, close all
%   parameter defination
global Tso Two Ms cp Mw t5 Ts5 Tsf
Vc = 0.50; %    cooler volume 0.50 ft^3
Vs = 0.15; %    volume of soda 0.15 ft^3 (12 cans)
Vw = Vc - Vs;
Two = 40; % initial temperature of water at 40 F
Tso = 80; % initial temperature of soda at 80 F
Ts5 = 68; % measured temperature of the soda after 5 min
Tsf = 55; % desired temperature of the soda is under 55 F
rho = 62.3; %   density of both water and soda 62.3 lbm/ft3 
cp = 1.0; % specific heat of both water and soda
Ms = rho*Vs;
Mw = rho*Vw;
to = 0; %   initial time is 0s
t5 = 5; %    measured temperature for sode at time is 5 min
tf = 40; %   set a final time of 40 min

%   using analytical solution to solve the hA value using fzero
% hA = 0.7321


hA1 = 1;
hA = fzero('Project1_P1_fzero',hA1)

%   using the analytically hA to get kw and ks
global kw ks
kw = hA/(Mw*cp);
ks = hA/(Ms*cp);

%   using the analytical solution to get the time to reach desired T of 55
%   tfe =19.9565

tfe = fzero('Project1_P1_fzero2',t5);
disp('the time from analytical solution to hit the desired Ts'), tfe

%   using analytical solution to get soda temperature vs time
Nt = 1000;
t = linspace(to,tf,Nt);
Ts = zeros(size(t));
Ts(1) = Tso;
Tw = zeros(size(t));
Tw(1) = Two;
for i = 2:Nt
    ti = t(i);
    Ts(i) = (Tso*kw+Two*ks)/(kw+ks)-(Two-Tso)/(kw+ks)*ks*exp(-(ks+kw)*ti);
    Tw(i) = (Tso*kw+Two*ks)/(kw+ks)+(Two-Tso)/(kw+ks)*kw*exp(-(ks+kw)*ti);
end


%   numerical slution to solve hA
%   setup constant matrices for state space formulation
%   looping HA to get numerical guess of HA
To = [Tso Two]';
tol = 1e-6;    
options = odeset('RelTol',tol);
j = 1;
n = 50000;
HA=0.0001;
error = 1;
tol2 = 1e-3;
while abs(error)>tol2 && j <= n
    KW = HA/(Mw*cp);
    KS = HA/(Ms*cp);
    KK = [-KS KS; KW -KW];
    [t1,T1] = ode23('Project1_P1_ode',[to t5],To,options,KK);
    e1 = T1(end,1)-Ts5;
    error = e1;
    if abs(error) > tol
        HA = HA+0.0001;
        j = j+1;
    end
end
disp('hA from numerical iterative loop is'),HA


[t2,T2] = ode23('Project1_P1_ode',[to tf],To,options,KK);

PP=spline(t2,T2(:,1));
global Tsf PP
tfd=fzero('tfe_P1',10)

disp('the time from numerical iterative loop to hit the desired Ts is'),tfd
%   plot analytical solution
plot(t,Ts,'r',t,Tw,'b',t2,T2(:,1),'k--*','LineWidth',2),
hold on
plot(t2,T2(:,2),'k--o')
title('Analytical Solution vs. Numerical Solution for Temperature of water and soda')
xlabel('time (min)'),ylabel('Temperature (deg F)'),grid
legend('Tsoda(exact)','Twater(exact)','Tsoda(ode23)','Twater(ode23)')





