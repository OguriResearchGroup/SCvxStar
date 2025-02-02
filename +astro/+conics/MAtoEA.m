function [EA, convflag] = MAtoEA(MA, e, tol)
    if nargin == 2
        tol = 1e-12;
    end        
    maxiter = 10;
    convflag = 0;

    x0 = MA; % initial guess
    for i = 1:maxiter
        f = x0-e*sin(x0)-MA;
        df = 1-e*cos(x0);
        x1 = x0 - f/df;
        if abs(f) < tol
            EA = x1;
            convflag = 1;
            return
        else
            x0 = x1;
        end
    end
    error('MA to EA conversion did not converge')
end