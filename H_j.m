%Enthalpy equation for the jth stage
%takes X, Feed - F and stage index j, and return jth stage balance
function [ H ] = H_j(X,j)
    c=5;
    n=19;
    B=62;
    reflux=9.5;
    H = zeros(c,1);
    L_sum = zeros(1,3);
    V_sum = zeros(1,3);
    h_Liq = zeros(1,3);
    h_Vap = zeros(1,3);
    if j==1
        H = X(c+2:2*c+1,1) - (reflux/(reflux+1))*X(1:c,1);
    elseif j==19
        for i=1:c
            H(i)= X(c+1+i,n)-B;
        end
    else
        for i=j-1:j+1
            L_sum(i) = sum(X(c+2:2*c+1, i));
            V_sum(i) = sum(X(1:c,i));
            h_Liq(i) = H_L(X(c+2:2*c+1, i)./L_sum(i) , X(c+1,i));
            h_Vap(i) = H_V(X(1:c,i)./V_sum(i) , X(c+1,i));
        end
        H = L_sum(j)*h_Liq(j) + V_sum(j)*h_Vap(j) -V_sum(j+1)*h_Vap(j+1) - L_sum(j-1)*h_Liq(j-1);
    end
end

