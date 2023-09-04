%
function Ts = fzero1(hA)
global Tso Two ms cp mw t5 Ts5 ks kw
ks = hA/ms*cp;
kw = hA/mw*cp;
Ts = ((Tso*kw+Two*ks)/(ks+kw))+(ks*(Tso-Two)*(exp(-(ks+kw)*t5))/(ks+kw))-Ts5;
end
%
% end of function