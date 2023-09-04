%   Defining the exact solution function of soda temperature
%
function Tse = fzero2(t)
global Tso Two kw ks Tsd
Tse = ((Tso*kw+Two*ks)/(ks+kw))+(ks*(Tso-Two)*(exp(-(ks+kw)*t))/(ks+kw))-Tsd;
end
%
% end of function