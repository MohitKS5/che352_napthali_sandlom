%% takes x and T and returns activity coefficients
function [gamma]=wilson(T,x)
    %given
    c=5;
    aij=[
    0.0             6.6429*10^2     8.34583*10^2    1.67946*10^3    1.70257*10^3
    -2.1501*10^2    0.0             -6.5210*10      4.949199*10^2   -7.3150*10
    -9.8420*10      1.6126*10^2     0.0             1.01700*10      7.91000*10
    2.1613*10^2     -1.6790*10^2    2.0054*10^2     0.0             1.4844*10^2
    -3.7225*10^2    -3.1310*10^2    -4.3141*10^2    -2.0850*10^2    0.0
    ];
    % from nist
    V=[40.73 74.05 79.84 89.41 80.67].*10^-3; %in m^3/mol
    
    R=1.987;
    gamma=zeros(1,5);
    A=ones(5,5);
    for i=1:c
        for j=1:c
            A(i,j)=(V(j)/V(i))*exp(aij(i,j)/(R*T));
        end
    end
    summation=A*x;
    for  i=1:5
        t = sum(x.*A(:,i)./summation);
        gamma(i)=exp(1-log(summation(i))-t);
    end
end

