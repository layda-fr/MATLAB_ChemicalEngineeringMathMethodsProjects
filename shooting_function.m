% Shooting Method
%
%%
function zp = shooting_function(x,z)
    global Tinf epsilon sigma Tsur P k Ac
%
    zp = zeros(2,1);
    hc = H_C(z(1));
    hr = H_R(z(1));
    zp(1)  = z(2);
    zp(2)  = ((hc*P)/(k*Ac))*(z(1)-Tinf)-((hr*P)/(k*Ac))*(z(1)-Tsur);  
%
    function hc = H_C(T)
        hc = 2.89*(0.6+0.624*(T-Tinf)^(1/6))^2;
    end
%
    function hr = H_R(T)
        hr = epsilon*sigma*(T+Tsur)*(T^2+Tsur^2); 
    end
end
%      