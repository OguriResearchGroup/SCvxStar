function TA = MAtoTA(MA,e,tol)
    import astro.conics.* % import astro package
    switch nargin
        case 2
            [EA, convflag] = MAtoEA(MA,e);
        case 3
            [EA, convflag] = MAtoEA(MA,e,tol);
        otherwise
            error('Number of input arguments must be 2 or 3.')
    end

    if convflag == 0
        TA = NaN;
        return
    end
    TA = EAtoTA(EA,e);
end