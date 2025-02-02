function cart_vec = mee2cart(x, MU)
    p = x(1);
    f = x(2);
    g = x(3);
    h = x(4);
    k = x(5);
    L = x(6);
    alpha = sqrt(h^2-k^2);
    s = sqrt(1+h^2+k^2);
    w = 1+f*cos(L)+g*sin(L);
    r = p/w;
    rvec = [r/s^2*(cos(L)+alpha^2*cos(L)+2*h*k*sin(L))
            r/s^2*(sin(L)-alpha^2*sin(L)+2*h*k*cos(L))
            2*r/s^2*(h*sin(L)-k*cos(L))];
    vvec = [-1/s^2*sqrt(MU/p)*(sin(L)+alpha^2*sin(L)-2*h*k*cos(L)+g-2*f*h*k+alpha^2*g)
            -1/s^2*sqrt(MU/p)*(-cos(L)+alpha^2*cos(L)+2*h*k*sin(L)-f+2*g*h*k+alpha^2*f)
            2/s^2*sqrt(MU/p)*(h*cos(L)+k*sin(L)+f*h+g*k)];
    cart_vec = [rvec; vvec];
end
