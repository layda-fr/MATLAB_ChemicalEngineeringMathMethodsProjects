% Simulation of a Simple Mechanical System
% Example for Numerical Integration of System of ODEs (IVPs)
%
% A simple mechanical system consisting of a mass, spring, and dashpot as well as a specified forcing function, f(t).  
% 
% A force balance on the system gives the following second order ODE:
%       my''(t) + cy'(t) + ky(t) = f(t)
%
%   For specificity, let
%       m =1,   k = 2,   c = 3,   and   f(t) = 10exp(-2t)
%   and for initial conditions, let
%       y(0) = 0   and   y'(0) = -5
%   all with appropriate units, of course.
%
%   This file also uses the following two files:
%       ivpnex1a.m - function file for ode23
%       ivpnex1b.m - function file for generating value of forcing function
%

      clear all, close all,  nfig = 0;
%
%   setup problem specific data
      m = 1;   k = 2;   c = 3;
      yo = 0;   ypo = -5;
%
%   setup constant matrices for state space formulation
      AA = [0 1; -k/m -c/m];   BB = [0 1/m]';
%
%   time domain of interest
      to = 0;   tf = 3;   te = linspace(to,tf,61);   
%
%   evaluate and plot forcing function over interval
      ff = ivpnex1b(te);
      nfig = nfig+1;   figure(nfig)
      plot(te,ff,'r-','LineWidth',2)
      title('IVPnex1:  Forcing Function Used in Simulation')
      xlabel('time'),ylabel('f(t)'),grid
%
%   solve 2x2 system of ODEs (Initial Value Prob -- IVP)
      zo = [yo ypo]';     
      tol = 0.001;    options = odeset('RelTol',tol);
      [t,z] = ode23('ivpnex1a',[to tf],zo,options,AA,BB);
%
%   exact solution (from lecture notes)
      te = te';    ze = zeros(length(te),2);
      ze(:,1) = -(10*te + 5).*exp(-2*te) + 5*exp(-te);
      ze(:,2) =  20*te.*exp(-2*te) - 5*exp(-te);
%
%   compare exact and numerical solution 
      nfig = nfig+1;   figure(nfig)
      subplot(2,1,1),plot(te,ze(:,1),'g-',t,z(:,1),'ro','LineWidth',2),grid on
      title('IVPnex1:  ODE23 (pts) vs Exact Solution (line) for Mechanical System')
      xlabel('time'),ylabel('position')
      legend('Exact','ODE23','Location','SouthEast')
      subplot(2,1,2),plot(te,ze(:,2),'g-',t,z(:,2),'ro','LineWidth',2),grid on
      title('IVPnex1:  ODE23 (pts) vs Exact Solution (line) for Mechanical System')
      xlabel('time'),ylabel('velocity')
      legend('Exact','ODE23','Location','SouthEast')

%

