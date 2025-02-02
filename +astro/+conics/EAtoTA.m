function TA = EAtoTA(EA,e)
    TA = 2*atan(sqrt((1+e)/(1-e)) * tan(EA/2));
end
