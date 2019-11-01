%% fugacity coefficient using Redlich Kwong Equation of state
% T in K and x[c]
function phi = fugacity(T,x)
    c=5;
    P = 101325;
    % NIST database critical T and P
    Pc = [81 48 47.5 48.9 53.2868]*10^5; %in Pa
    Tc = [513 508 510 562 537]; %in Kelvin

    % parameters of pure component Redlich Kwong Equation of State
    ai = zeros(5,1);
    bi = zeros(5,1); 
    Ai = zeros(5,1);
    Bi = zeros(5,1);

    %gas constant
    R = 8.314;

    for i = 1:c
        ai(i) = (0.42478 * R * R * (Tc(i)^2.5))/Pc(i);
        Ai(i) = (ai(i)^0.5)/(R * (T^1.25));
        bi(i) = 0.08664 * R * Tc(i) / Pc(i);
        Bi(i) = bi(i)/(R * T);
    end
    
    % parameters of the mixture
    A = sum(x.*Ai);
    B = sum(x.*Bi);

    % solve for Z
    polyn = [1,-1,((A^2)*P-B*P*(1 + B*P)),(-1*(A*P)^2*B)];
    r = roots(polyn);
    % extract the real positive(so max) value
    r =r(imag(r)<10^-10);
    z=max(r);
    % calc fugacities
    phi = zeros(1,c);
    for j=1:c
        phi(j) = exp(Bi(j)*(z-1)/B - log(z - B*P) - A^2/B*(2*(Ai(j)/A) - Bi(j)/B)*log(1 + (B*P/z)));
    end
end
