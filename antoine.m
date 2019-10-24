 %% takes (T,P) and returns (Psat, k)
function [Psat,k] = antoine(T)
    c=5;
    P=1.01325;
    
    % from nist database
    A=[5.20409 4.42448 4.20364 4.72583 4.20772];
    B=[1581.341 1312.253 1164.426 1660.652 1233.129];
    C=[-33.50 -32.445 -52.69 -1.461 -40.953];
    
    Psat=zeros(1,c);
    k=zeros(1,c);
    for i=1:c
        Psat(i)=10^(A(i)-B(i)/(T+C(i)));
        k(i)=Psat(i)/P;
    end
end
