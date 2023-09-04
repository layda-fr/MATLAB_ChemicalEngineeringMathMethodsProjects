%%
% FD Method
%
%%
%   getting started 
clear all, close all; nfig=0; clc
%
%% Problem Variables
%   
global Tb Tinf Tsur k Ac hc hr L D P epsilon sigma
%
%   Initially given variables
%                   % UNIT                 DESCRIPTION
%                                   
Tb = 100;           % oC                   The fixed base temperature
L = 0.1;            % m                    The length of the fin
epsilon = 0.35;     
Tinf = 25;          % oC                   The temperature at infinity
D = 0.01;           % m                   
sigma = 5.67e-8;    % W/m^2-oK^4         
k = 14;             % W/m-oK               
Tsur = 25;          % oC                  The temperature of the surface
%
%  Other variables that have been defined
%
P=D*pi;                             % The perimeter of the circular cross section of the pin-shaped fin in meters
Ac=(D^2)*pi/4;                      % The area of the circular cross section of the pin-shaped fin in meters
Acl=Ac*L;                           % The surface area of the circular cross section of the pin-shaped fin in meters
xo=0;   
xf=L; 
To=Tb;
%
%   user defined number of unknowns in problem
N = 25;
%
% compute interval size and discrete x vector
x = linspace(xo,xf,N);
dx = (xf-xo)/(N-1); 
x = x';
%
yxo = Tb;
% specify yold
yold = yxo*ones(size(x));
%
% start iteration loop
itmax = 200; 
it = 0; 
tol = 1e-6; 
emax = 1;
while emax > tol && it <= itmax
it = it+1; 
%
% compute solution profile 
[a,b] = FD_2(yold,dx);                   % determine the coefficient matrices
ynew = a\b;                                 % find the new solution vector
%
% check convergence
emax = max(abs((ynew-yold)./ynew));
%
% define guess for next iteration
yold = ynew;
end
%
% add boundary points to solution for plotting
yn = [yold]; xn = [x];
%
%
qconv = (H_C(yold(1:end-1,1)).*(yold(1:end-1,1)-Tinf).*(L*P+Ac*L))*dx;   
qconv = sum(qconv);
qconv = qconv + H_C(yold(end,1))*Ac*(yold(end,1)-Tinf);
disp(['qconv = ', num2str(qconv)]);
%
qrad = (H_R(yold(1:end-1,1)).*(yold(1:end-1,1)-Tsur).*(L*P+Ac*L))*dx;
qrad = sum(qrad);
qrad = qrad + H_R(yold(end,1))*Ac*(yold(end,1)-Tsur);
disp(['qrad = ', num2str(qrad)]);
%
qtot = qconv + qrad;
disp(['qtot = ', num2str(qtot)]);
%
% plot final results
nfig = nfig+1; figure(nfig)
plot(xn,yn,'r-','LineWidth',2),hold on
title('FD Method Solution') 
xlabel('X values'),ylabel('Temperature values'),grid
legend('FD','Location','SouthEast')
%
    function hc = H_C(T)
        hc = 2.89.*(0.6+0.624.*(T-25).^(1/6)).^2;
    end
    function hr = H_R(T)
        hr = 0.35.*5.67e-8.*(T+25).*(T.^2+25^2); 
    end
    function E = E_1(T)
        E = H_C(T).*P./(k.*Ac);
    end        
    function E = E_2(T)
        E = H_R(T).*P./(k.*Ac);          
    end
    function R = R_1(T)
        R = H_C(T)./(k);
    end
    function R = R_2(T)
        R = H_R(T)./(k);          
    end
%
% end
