classdef Constants

    properties (Constant)

        %% Time Constants:
        HOUR2SEC = 3600;
        DAY2SEC = 86400;
        YEAR2DAY = 365.25;
        YEAR2SEC = 86400 * 365.25;
        SEC2HOUR = 1/3600;
        SEC2DAY = 1/86400;
        SEC2YEAR = 1/(86400 * 365.25);

        %% Celestial Body Radii:
        % Data from JPL/NAIF's PCK00010.TPC Planetary Constants Kernel.

        SUN_radii_KM        = 696000.00000;
        MERCURY_radii_KM    = 2439.70000; 
        VENUS_radii_KM      = 6051.80000; 
        EARTH_radii_KM      = 6378.13660; 
        MOON_radii_KM       = 1737.40000; 

        %% Celestial Body Gravitational Parameters:
        % Data from JPL/NAIF's DE440.BSP Kernel, summarized here:
        % https://ssd.jpl.nasa.gov/astro_par.html

        SUN_GM_KM3PS2        = 132712440041.939400;
        MERCURY_GM_KM3PS2     = 22031.868551;
        VENUS_GM_KM3PS2       = 324858.592000;
        EARTH_GM_KM3PS2       = 398600.435507;
        MOON_GM_KM3PS2        = 4902.800118;
        MARS_SYS_GM_KM3PS2    = 42828.375816;
        JUPITER_SYS_GM_KM3PS2 = 126712764.100000;
        SATURN_SYS_GM_KM3PS2  = 37940584.841800;
        URANUS_SYS_GM_KM3PS2  = 5794556.400000;
        NEPTUNE_SYS_GM_KM2PS3 = 6836527.100580;
        PLUTO_SYS_GM_KM2PS3   = 975.500000;

        AU = 1.49597870700E+8; % [km]

    end

end