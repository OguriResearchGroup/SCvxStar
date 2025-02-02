% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % EarthMars.m returns a struct with parameters for the Earth-Mars transfer problem
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = EarthMars()

    import astro.Constants;
    import astro.ScaleFactor;
    import astro.TwoBody;

    p = minfuelparams();

    p.CntrlType            = 'LowThrust';
    p.Nseg                 = 30;
    p.tof_days             = 500;
    p.et0                  = cspice_str2et('2024 Aug 11 12:00:00');
    p.Spacecraft.Tmax_dim  = 0.5 * 1E-03; % [kN]
    p.Spacecraft.mass0_dim = 1000; % [kg]

    % ScaleFactor and DynamicalSystem
    LU = Constants.AU;
    TU = sqrt(LU^3/Constants.SUN_GM_KM3PS2);
    MU = p.Spacecraft.mass0_dim;
    SF = ScaleFactor(LU,TU,MU);
    p.DS = TwoBody(1,SF);

    p.tof = p.tof_days * SF.DAYS2TU;
    p.t_his = linspace(0,p.tof, p.Nseg+1);

    tof_dim = p.tof * SF.t;
    p.etf = p.et0 + tof_dim;

    X_SF = [SF.l * ones(3,1); SF.v * ones(3,1)];

    SPICE_FRAME    = 'ECLIPJ2000';
    SPICE_ABCORR   = 'none';
    SPICE_OBSERVER = 'Sun';

    p.x_init = cspice_spkezr('EARTH', p.et0, SPICE_FRAME, SPICE_ABCORR, SPICE_OBSERVER)./X_SF;
    p.x_fin = cspice_spkezr('MARS', p.etf, SPICE_FRAME, SPICE_ABCORR, SPICE_OBSERVER)./X_SF;

    p.u_max = (p.Spacecraft.Tmax_dim/SF.force) / (p.Spacecraft.mass0_dim / SF.mass);
end