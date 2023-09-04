% Simulation of the Dynamics of a Semi-Batch Chemical Reactor
% A Typical Initial Value Problem                         
%                                                                               
% This project involves the modeling/solution of the dynamics of a semi-batch chemical reactor.  
% A semi-batch reactor is one that has flow in, but no flow out.
% The specific system of interest has three components and the reactor is fed with constant flow and feed concentration of species A.  Components B and C  are formed from component A.
%
% The defining eqns for the number of moles of each component are:
%      dNA/dt = Qo*CAo - (k1/VR)*NA^2 + k2*NB
%      dNB/dt = (k1/VR)*NA^2 - k2*NB - (k3/VR)*NB^2
%      dNC/dt = (k3/VR)*NB^2
% with
%      VR = Vo + Qo*t
% and initial conditions
%      NAo = NBo = NCo = 0    (all components initially at zero concentration)
%
%   These equations will be solved for a specified set of data.  The goal is
%   to determine the concentrations versus time, Ci = Ni/VR, with particular
%   focus on the time when the second component (B) reaches its maximum
%   concentration.
%
      clear all, close all, nfig = 0;
%
%   generic and plant specific data
      Qo = 10.0;      % inlet flow rate (liters/sec)
      Vo = 50.0;      % initial volume of semi-batch reactor (liters)
      CAo = 1.0;      % concentration of component A in feed (gmole/liter)
      K1 = 0.1;  K2 = 0.05;  K3 = 0.05;  % rate constants 
%
%   now simulate system for several seconds (past peak in concentration in B)
      to = 0;    tf = 100;    % time interval 
      no = [0 0 0]';          % initial conditions 
      tol = 0.001;    options = odeset('RelTol',tol);
      [t,n] = ode23('semibra',[to tf],no,options,Qo,Vo,CAo,K1,K2,K3); 
%
%   compute concentrations
      VR = Vo + Qo*t;    C = zeros(size(n));
      for i = 1:3    C(:,i) = n(:,i)./VR;    end
%
%   plot primary results 
      nfig = nfig+1;  figure(nfig)
      plot(t,C(:,1),'b-',t,C(:,2),'r-.',t,C(:,3),'g--','LineWidth',2),grid
      title('SemiBR:  Component Concentrations (gmoles/liter)')                                          
      xlabel('Time (seconds)'),ylabel('Concentrations (gmoles/liter)')
      legend('Component A','Component B','Component C')                   
%                                                                              

