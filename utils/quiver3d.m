function quiver3d(x, u, varargin)
    quiver3(x(1,:), x(2,:), x(3,:), u(1,:), u(2,:), u(3,:), varargin{:});
end