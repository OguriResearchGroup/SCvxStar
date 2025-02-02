function coe = cart2coe(rv,mu,anomaly)
    arguments
        rv (6,1) double
        mu (1,1) double
        anomaly string = 'TA'
    end

    eps = 1.e-10;
    R = rv(1:3);
    V = rv(4:6);
    r = norm(R);
    v = norm(V);
    vr = dot(R,V)/r;
    H = cross(R,V);
    h = norm(H);
    %...Equation 4.7:
    incl = acos(H(3)/h); % There is no quadrant ambiguity here! (Curtis)
    %...Equation 4.8:
    N = cross([0 0 1],H);
    n = norm(N);
    %...Equation 4.9:
    if n == 0
        W = 0;
        warning('This orbit has inclination close to 0 or pi. May cause singularities.')
    else
        W = atan2(N(2),N(1)); %(From yuri's notes)   
    end
    %...Equation 4.10:
    E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
    e = norm(E);
    %...Equation 4.12 (incorporating the case e = 0):
    if n ~= 0
        if e > eps
            w = acos(dot(N,E)/n/e);
            if E(3) < 0
                w = 2*pi - w;
            end
        else
            w = 0;
        end
    else
        w = 0;
    end
    %...Equation 4.13a (incorporating the case e = 0):
    if e > eps
    %     TA = atan2(h*vr, h^2/r-mu); %Eq 3.29 of Battin)
    % E
    % R
    % quiver(0,0,E(1),E(2))
    % hold on
    % quiver(0,0,R(1),R(2))
        TA = acos(dot(E,R)/e/r);
        if vr < 0
            TA = 2*pi - TA;
        end
    else
        warning('This orbit has eccentricity close to 0. May cause singularities.')
        cp = cross(N,R);
        if cp(3) >= 0
            TA = acos(dot(N,R)/n/r);
        else
            TA = 2*pi - acos(dot(N,R)/n/r);
        end
    end
    %...Equation 4.62 (a < 0 for a hyperbola):
    a = h^2/mu/(1 - e^2);

    coe = [a e incl W w TA]';
    switch anomaly
        case 'MA'
            coe(6) = TAtoMA(coe(6),e);
        case 'EA'
            coe(6) = TAtoEA(coe(6),e);
    end


end