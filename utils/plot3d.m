function varargout = plot3d(x, varargin)
    [varargout{1:nargout}] = plot3(x(1,:), x(2,:), x(3,:), varargin{:});
    axis equal
end