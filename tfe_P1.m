function fh = tfe_P1(t)
global Tsf PP
fh = ppval(PP,t)-Tsf;