%%
%
function ti = fzero3(t)
global Tsd PP
ti = ppval(PP,t)-Tsd;
end
%
% end of function