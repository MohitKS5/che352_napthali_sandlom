%% takes F,k[c] and z[c] array and returns V,L,x[c],y[c]
function [L,V,x,y] = rachford_rice(F,k,z)
    c = length(k);
    f=@(phi) 0;
    for i=1:c
        f = @(phi) f(phi)+z(i)*(k(i)-1)/(1+phi*(k(i)-1));
    end
    p = newton_raphson_1d(f,0.5);
    
    V = p*F;
    L = (1-p)*F;
    % calculate compositions
    x=zeros(c,1);
    y=zeros(c,1);
    for i=1:c
        x(i) = F*z(i)/(L+k(i)*V);
        y(i) = k(i)*x(i);
    end
end
