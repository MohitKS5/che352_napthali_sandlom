function [ dep ] = dep_L(x, T_stage )
    c=5;
    P = 101325;
    % NIST database critical T and P
    Pc = [81 48 47.5 48.9 53.2868]*10^5; %in Pa
    Tc = [513 508 510 562 537]; %in Kelvin
    Vm = [40.73; 74.05; 79.84;89.41;80.67].*10^-6; %in m^3/mol
    % parameters of pure component Redlich Kwong Equation of State
    ai = zeros(5,1);
    bi = zeros(5,1); 
    %gas constant
    R = 8.314;
    for i = 1:c
        ai(i) = (0.42478 * R * R * (Tc(i)^2.5))/Pc(i);
        bi(i) = 0.08664 * R * Tc(i) / Pc(i);
    end
    V = 0;
    for i=1:c
        V = V + x(i)*Vm(i);
    end;
    z = P*V/(R*T_stage);
    a_m = 0;
    b_m = 0;
    for i=1:c
        for j=1:c
            a_m + a_m + x(i)*x(j)*((ai(i)*ai(k))^0.5);
        end
        b_m = b_m + x(i)*bi(i);
    end
    dep = R*T_stage*(-1 + z -(3*a_m*log(1+b_m/z)/(2*b_m)));
end

