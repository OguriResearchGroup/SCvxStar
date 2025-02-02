function EA = TAtoEA(TA, e)
    EA = 2*atan(sqrt((1-e)/(1+e)) * tan(TA/2));
end