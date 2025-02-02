function mee = coe2mee(coe)
    a = coe(1);
    e = coe(2);
    i = coe(3);
    W = coe(4);
    w = coe(5);
    TA = coe(6);
    p = a*(1-e^2);
    f = e*cos(w+W);
    g = e*sin(w+W);
    h = tan(i/2)*cos(W);
    k = tan(i/2)*sin(W);
    L = W+w+TA;
    mee = [p f g h k L]';
end