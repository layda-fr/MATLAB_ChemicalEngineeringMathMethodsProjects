%   given differential eqn.   z' = Az + g   where g = B*u
%
%   uses function file ivpnex1b.m to evaluate forcing function at time t
%
      function zp = ivpnex1a(t,z,flag,AA,BB)
      u=ivpnex1b(t);
      g = BB*u;
      zp = AA*z+g;
      end
%
%   end of function
