%   The defining eqns are:
%      dNA/dt = Qo*CAo - (k1/VR)*NA^2 + k2*NB
%      dNB/dt = (k1/VR)*NA^2 - k2*NB - (k3/VR)*NB^2
%      dNC/dt = (k3/VR)*NB^2
%   with
%      VR = Vo + Qo*t
%
       function np = odefile(t,n,flag,Qo,Vo,CAo,K1,K2,K3)
       np = zeros(length(n),1);
       VR = Vo + Qo*t;
       np(1) = Qo*CAo - (K1/VR)*n(1)^2 + K2*n(2);
       np(2) = (K1/VR)*n(1)^2 - K2*n(2) - (K3/VR)*n(2)^2;
       np(3) = (K3/VR)*n(2)^2;
%
%   end of function

