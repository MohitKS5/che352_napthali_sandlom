function root = newton_raphson_1d(f,x0)
    error = 10^(-4);
    syms x;
    d= matlabFunction( diff(f(x)));
    if nargin(d) == 0 
        d = @(x) d();
    end
    
    x1=x0;
    for i=1:100
        x2=x1-(f(x1))/(d(x1));
        err=abs((x2-x1)/x1);
        if err<error
            break
        end
        x1=x2;
    end
    root=x1;
end