% Shooting Method
%
%%
%   getting started
clear all, close all; nfig=0; clc
%
%% Problem Variables
global Tb Tinf Tsur k Ac L D P epsilon sigma
%   Initially given variables
%                   % UNIT          %       DESCRIPTION
%                                   %       
Tb = 100;           % oC            %       The fixed base temperature
L = 0.1;            %  m            %       The length of the fin
epsilon = 0.35;     
Tinf = 25;          % oC            %       The temperature at infinity
D = 0.01;           % m             %       
sigma = 5.67e-8;    % W/m^2-oK^4    %      
k = 14;             % W/m-oK        %       
Tsur = 25;          % oC            %      The temperature of the surface
%
%  Other variables that have been defined
%
P=D*pi;                             % The perimeter of the circular cross section of the pin-shaped fin in meters
Ac=D^2*pi/4;                        % The area of the circular cross section of the pin-shaped fin in meters
Acl=Ac*L;                           % The surface area of the circular cross section of the pin-shaped fin in meters
xo=0;   
xf=L; 
To=Tb;
%
%
%% Shooting Method
    %
    %
    tol     = 1e-6;  
    options = odeset('RelTol',tol);    % set tolerance for ode45 
    [xs,zs] = mybvpsh('shooting_function',[xo xf],options);   
    dx = L / length(zs(:,1));
%
qconv = (H_C(zs(1:end-1,1)).*(zs(1:end-1,1)-Tinf).*(L*P+Ac*L))*dx;   
qconv = sum(qconv);
qconv = qconv + H_C(zs(end,1))*Ac*(zs(end,1)-Tinf);
disp(['qconv = ', num2str(qconv)]);
%
qrad = (H_R(zs(1:end-1,1)).*(zs(1:end-1,1)-Tsur).*(L*P+Ac*L))*dx;
qrad = sum(qrad);
qrad = qrad + H_R(zs(end,1))*Ac*(zs(end,1)-Tsur);
disp(['qrad = ', num2str(qrad)]);
%
qtot = qconv + qrad;
disp(['qtot = ', num2str(qtot)]);
%
%   plot final function values 
      plot(xs,zs(:,1),'LineWidth',2),grid
      legend('Shooting Method','Location','SouthEast')
      title('Shooting Method for ODEs with Automated Iteration')                  
      xlabel('X values'),ylabel('Temperature values')                        
%   
    function hc = H_C(T)
        hc = 2.89.*(0.6+0.624.*(T-25).^(1/6)).^2;
    end
%
    function hr = H_R(T)
        hr = 0.35.*(5.67e-8).*(T+25).*(T.^2+25^2); 
    end
% end 

