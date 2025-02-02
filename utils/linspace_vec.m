function x_interp = linspace_vec(x1, x2, n)
    assert(length(x1) == length(x2), 'x1 and x2 must have the same length');
    x_interp = NaN(length(x1), n);
    for i = 1:length(x1)
        x_interp(i,:) = linspace(x1(i), x2(i), n);
    end
end