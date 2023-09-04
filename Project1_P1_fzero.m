%   define the fzero function for the exact solution of soda temperature
function y=Project1_P1_fzero(hA)
global Tso Two Ms cp Mw t5 Ts5
kw = hA/(Mw*cp);
ks = hA/(Ms*cp);
y=(Tso*kw+Two*ks)/(kw+ks)+(Tso*ks-ks*Two)/(kw+ks)*exp(-(ks+kw)*t5)-Ts5;