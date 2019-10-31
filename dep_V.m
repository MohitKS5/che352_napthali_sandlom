function [ dep ] = dep_V( y, T_stage )
    c=5;
    P = 101325;
    % NIST database critical T and P
    Pc = [81 48 47.5 48.9 53.2868]*10^5; %in Pa
    Tc = [513 508 510 562 537]; %in Kelvin
    % parameters of pure component Redlich Kwong Equation of State
    ai = zeros(5,1);
    bi = zeros(5,1); 
    %gas constant
    R = 8.314;
    for i = 1:c
        ai(i) = (0.42478 * R * R * (Tc(i)^2.5))/Pc(i);
        bi(i) = 0.08664 * R * Tc(i) / Pc(i);
    end
    a_m = 0;
    b_m = 0;
    for i=1:c
        for j=1:c
            a_m + a_m + x(i)*x(j)*((ai(i)*ai(k))^0.5);
        end
        b_m = b_m + x(i)*bi(i);
    end
    t1 = -1*R*T_stage;
    t2 = a_m/(T_stage^0.5) - R*T*b_m - P*(b_m^2);
    t3 = -1*(a_m*b_m)/(T_stage^0.5);
    e = [P t1 t2 t3];
    r = roots(e);
    r = r(imag(r) == 0);
    V = max(r);
    z = (P*V)/(R*T);
    dep = R*T_stage*(-1 + z -(3*a_m*log(1+b_m/z)/(2*b_m)));
end

