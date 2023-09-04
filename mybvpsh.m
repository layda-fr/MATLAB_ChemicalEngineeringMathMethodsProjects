% Shooting Method
%
%%
function [t,z]  = mybvpsh(func,tspan,options)

global Tb epsilon sigma Tsur Tinf k
 alfo = 1;
    %
    %   extract BC coefficients (for convenience)
      
      zo=[Tb;alfo];
    %
    %   set some other iterative parameters
      err = 1e10;   tol = 1e-8;   icnt = 1;   icntmax = 25;
    %
    %   start iteration  
      while abs(err) > tol  &&  icnt <= icntmax
                 zo     = [Tb; alfo];
                [t,z]   = ode23(func,tspan,zo,options);
                e1      = z(end,2)-BV(z(end,1));
                err     = e1;
                if abs(err) > tol
%                    Our count exceeded is only checked if error is above
%                    tol & gives a message if it happens. 
                    if icnt>icntmax
                        disp('itm surpassed');
                        break
                    end
                    alfp    = 1.01*alfo;
                    zo      = [Tb; alfp];
                   [t,z]    = ode23(func,tspan,zo,options);
                    e2      = z(end,2)-BV(z(end,1));
                end
                deda = (e2-e1)/(0.01*alfo);    
                alfn = alfo-e1/deda;
                icnt = icnt+1;  
                alfo = alfn;  
      end   % end of iteration loop
%
      if icnt >= icntmax
        fprintf(1,'  WARNING -- Hit maximum iteration limit!!! \n');
      end
%  
    function hc = H_C(T)
        hc = 2.89*(0.6+0.624*(T-Tinf)^(1/6))^2;
    end
%
    function hr = H_R(T)
        hr = epsilon*sigma*(T+Tsur)*(T^2+Tsur^2); 
    end
%
    function BV = BV(T)
        term1   = H_C(T)/(-k)*(T-Tinf);
        term2   = H_R(T)/(-k)*(T-Tsur);            
        BV      = term1 + term2;
    end
end
