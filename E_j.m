% material eq for the jth tray
% takes X , Feed and stage index j and returns jth stage balance
function E =E_j(X,j)
    %% given
    c=5;
    P = 1.01325;
    
    %% get mole fractions in gas and liquid
    v = X(1:c,j);
    l = X(c+2:2*c+1,j);
    y=v./sum(v);
    x=l./sum(l);
    %% Equations E = 
    fug = fugacity(X(c+1,j),x);
    activity = wilson(X(c+1,j),x);
    [psat,~] = antoine(X(c+1,j));
    E = P*(fug.*y') - activity.*psat.*x';
end
