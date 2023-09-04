%   define the fzero function for the exact solution of soda temperature
function te=Project1_P1_fzero2(t)
global Tso Two kw ks Tsf
te=(Tso*kw+Two*ks)/(kw+ks)+(Tso*ks-ks*Two)/(kw+ks)*exp(-(ks+kw)*t)-Tsf;