function MA = TAtoMA(TA,e)
    import astro.conics.*; % import astro package
    EA = TAtoEA(TA,e);
    MA = EAtoMA(EA,e);
end