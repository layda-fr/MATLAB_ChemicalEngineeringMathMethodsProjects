%   PROJECT 1
%
% The project of "Keeping It Cool" focuses on a common occurrence of cooling a couple of 
% six-packs of soda by putting the twelve cans in a cooler loaded with chilly water from a 
% mountain range stream through a well-insulated cooler. In our situation, there will be no
% energy loss to the environment because the sole energy transfer will be between the soda and
% the water. The main purpose of this project is to estimate how long it will take to chill the 
% soda to a specifically desired temperature by exchanging energy with water, based on the 
% above situation. This aim is achieved by two different unforced linear time-invariant (LTI) 
% initial value problem (IVP)- based approaches: analytical solution via hand-written equations 
% by using eigenvalues, eigenvectors, and IVP; and numerical solution via Matlab by using a
% built-in ODE solver as ode45. 
%
%%
%   getting started
clear all, close all;
%
%% Problem Variables
%   
global Tso Two ms cp mw t5 Ts5 Tsd
%
%   Initially given variables
%                   % UNIT          %       DESCRIPTION
%                                   %
Vc = 0.50;          % ft3           %       volume of cooler
Vs = 0.15;          % ft3           %       volume of soda (12 cans)
Vw = Vc-Vs;         % ft3           %       volume of water
Two = 40;           % oF            %       initial temperature of water
Tso = 80;           % oF            %       initial temperature of soda
Ts5 = 68;           % oF            %       measured temperature of the soda after 5 minutes
Tsd = 55;           % oF            %       desired temperature
rho = 62.3;         % lbm/ft3       %       density of both water and soda
ms = rho*Vs;        % lbm           %       mass of both soda
mw = rho*Vw;        % lbm           %       mass of both water
cp = 1.0;           % Btu/lbm-oF    %       specific heat of both water and soda
%
%  Other variables that have been defined
%
to = 0;
t5 = 5;
tf = 60;
Na = 1000;
%
%% Analytical solution part
%% Evaluating hA by fzero function using the calculation of hA from exact solution
%
hA1 = 0.7266;
hA = fzero('fzero1',hA1);
%
%% Finding kw and ks by hA
% 
global kw ks
%
kw = hA/(mw*cp);
ks = hA/(ms*cp);
%
%% Displaying the found hA, kw, and ks values by exact solution
%
disp('hA found by exact solution'),hA
disp('kw found by hA and exact solution'),kw
disp('ks found by hA and exact solution'),ks
% 
%% Finding and displaying the time to reach desired temperature of soda
%
t = fzero('fzero2',t5);
disp('Time to reach desired temperature of soda'), t
%
%% Evaluating hA with Ts5 by using loop
%
t = linspace(to,tf,Na);
Twe = zeros(size(t));
Twe(1) = Two;
Tse = zeros(size(t));
Tse(1) = Tso;

for i = 2:Na
    te = t(i);
    Tse(i) = ((Tso*kw+Two*ks)/(ks+kw))-(ks*(Two-Tso)*(exp(-(ks+kw)*te))/(ks+kw));
    Twe(i) = ((Tso*kw+Two*ks)/(ks+kw))+(kw*(Two-Tso)*(exp(-(ks+kw)*te))/(ks+kw));
end
%
%% Numerical solution part
%
% We can calculate hA by given data point numerically. However, as time changes, hA would be
% changed. A numerical iteration to keep hA close to the numerical value found by analytical
% solution should be in place. Based on that, a new hA defined by initial HA. Besides, two 
% tolerance has been determined as tol1 = 1e-6 for ODE solver and tol2 = 1e-3 as threshold of
% while loop for evaluating original ODE. At the end of each loop, HA increased with 1e-4 and
% the next loop checked with the new HA value. When the desired HA value has been reached,
% the loop ends and displays the iterative result.
%
%% ODE Solver Setup     
%
% Setting iterative parameters
error = 1;  tol1 = 1e-6;  tol2 = 1e-3;  
icnt = 1;  icntmax = 1e6;
options = odeset('RelTol',tol1);
To = [Tso Two]';
HA = 1e-4;
%
% Starting iteration  
while abs(error) > tol2  &&  icnt <= icntmax
    KS = HA/ms*cp;
    KW = HA/mw*cp;
    A = [-KS KS;
         KW -KW];
    [t1,T1] = ode45('ode1',[to t5],To,options,A); 
    T_array = T1(end,1) - Ts5;
    error = T_array;
    if abs(error) > tol1
       HA = HA + 1e-4;
       icnt = icnt+1;
    end
end
%
disp('hA by numerical iteration'),HA
%
% Based on the iterative procedure, hA value is converged to the result found by analytical
% solution. 
%
%% Evaluation of the numerical solution with ODE45
%
[t2,T2] = ode45('ode1',[to tf],To,options,A);
%
%% Evaluating and displaying the time to reach desired temperature of soda
%
global Tsd PP
PP = spline(t2,T2(:,1));

ti = fzero('fzero3',19.9565);
disp('Time to reach desired temperature of soda by iteration' ), ti

%% Plot Information
%
nfig = 0;
nfig = nfig+1; figure(nfig)

subplot(2,2,1),plot(t2,T2(:,1),'b--s','LineWidth',2),grid on
title('Tsoda(t) Profile')
xlabel('Time (minutes)'),ylabel('Temperature (oF)')
legend('Numerical','Location','East')

subplot(2,2,2),plot(t,Tse,'r--o','LineWidth',2),grid on
title('Tsoda(t) Profile')
xlabel('Time (minutes)'),ylabel('Temperature (oF)')
legend('Analytical','Location','East')

subplot(2,2,[3,4]); plot(t2,T2(:,1),'b--s',t,Tse,'r--o','LineWidth',2),grid on
title('Tsoda(t) Profile for Both Numerical and Analytical Solutions' )
xlabel('Time (minutes)'),ylabel('Temperature (oF)')
legend('Numerical','Analytical','Location','East')


nfig = nfig+1; figure(nfig)

subplot(2,2,1),plot(t2,T2(:,2),'b--s','LineWidth',2),grid on
title('Twater(t) Profile')
xlabel('Time (minutes)'),ylabel('Temperature (oF)')
legend('Numerical','Location','East')

subplot(2,2,2),plot(t,Twe,'r--o','LineWidth',2),grid on
title('Twater(t) Profile')
xlabel('Time (minutes)'),ylabel('Temperature (oF)')
legend('Analytical','Location','East')

subplot(2,2,[3,4]); plot(t2,T2(:,2),'b--s',t,Twe,'r--o','LineWidth',2),grid on
title('Twater(t) Profile for Both Numerical and Analytical Solutions' )
xlabel('Time (minutes)'),ylabel('Temperature (oF)')
legend('Numerical','Analytical','Location','East')


nfig = nfig+1; figure(nfig)

subplot(2,2,[1,2]); plot(t2,T2(:,1),'b--s',t2,T2(:,2),'r--o','LineWidth',2),grid on
title('Overlapped Numerical Profiles of Tsoda(t) and Twater(t)' )
xlabel('Time (minutes)'),ylabel('Temperature (oF)')
legend('Numerical for Ts(t)','Numerical for Tw(t)','Location','SouthEast')

subplot(2,2,[3,4]); plot(t,Tse,'b--s',t,Twe,'r--o','LineWidth',2),grid on
title('Overlapped Analytical Profiles of Ts(t) and Tw(t)' )
xlabel('Time (minutes)'),ylabel('Temperature (oF)')
legend('Analytical for Tsoda(t)','Analytical for Twater(t)','Location','SouthEast')
%    
% end 
%
% The results of this project indicated that the numerical results of the computational
% method and the analytical results agree and converge with each other. The plot of the
% final converged Ts(t) and Tw(t)  profiles and the time it takes to cool the soda to 
% the desired temperature of 55 oF have been obtained.
%