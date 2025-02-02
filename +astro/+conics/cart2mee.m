function mee = cart2mee(varargin)
% Converts state vector to Modified Equinoctial Elements.
% Has an intermediate conversion to the Classical elements, which might be
% a problem
% Not used in Problem Set 1
    import astro.conics.*; % import astro package
    coe = cart2coe(varargin{:});
    mee = coe2mee(coe);
end