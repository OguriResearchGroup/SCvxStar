function labels3d(units)
    if nargin < 1
        units = 'ND';
    end

    xlabel(sprintf('X [%s]', units))
    ylabel(sprintf('Y [%s]', units))
    zlabel(sprintf('Z [%s]', units))
end