function coe = mee2coe(x)
    p = x(1);
    f = x(2);
    g = x(3);
    h = x(4);
    k = x(5);
    L = x(6);
    a = p/(1-f^2-g^2);
    e = sqrt(f^2+g^2);
    i = atan2(2*sqrt(h^2+k^2), 1-h^2-k^2);
    w = atan2(g*h-f*k, f*h+g*k);
    W = atan2(k, h);
    TA = L - atan2(g, f);
    coe = [a e i W w TA]';
end
