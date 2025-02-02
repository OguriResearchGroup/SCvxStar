function stairsZOH(t,x,varargin)
    assert(length(t)==length(x)+1)
    stairs(t, [vec(x); x(end)], varargin{:})
end