% FD Method
%
%%
function [A,B] = FD_2(T,h)
    % local variables 
    N = length(T);
    A = zeros(N,N);     %  initialize A to save memory reallocation. 
    B = zeros(N,1);     %  initialize B to save memory reallocation. 
                        %  Note we are setting up for the form
                        %  A(N,N)Y(N,1) = B(N,1)
%                       
global Tinf epsilon sigma Tsur P k Ac Tb
%
       E1=E_1(T);
       E2=E_2(T);
       R1=R_1(T);
       R2=R_2(T);
 %      
    % left boundary                                                             
      B(1)      = Tb;
      A(1,1)    = 1;                                                     
      A(1,2)    = 0;  
    % Central Points  
    for i = 2:N-1                                                          
      A(i,i-1)  = 1;    
      A(i,i)    = -(2+(E1(i)+E2(i))*(h^2));                                                     
      A(i,i+1)  = 1;     
      B(i)      = -(E1(i)*Tinf-E2(i)*Tsur)*(h^2);   
    end                                                                       
    %   right boundary                                                                  
    A(N,N-1)    = -1;                                          
    A(N,N)      = (1+h*(R1(N)-R2(N)));                                              
    B(N)            = (R2(N)*Tsur-R1(N)*Tinf)*(h);
%
    function hc = H_C(T)
        hc = 2.89.*(0.6+0.624.*(T-Tinf).^(1/6)).^2;
    end
    function hr = H_R(T)
        hr = epsilon.*sigma.*(T+Tsur).*(T.^2+Tsur^2); 
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
end 
%