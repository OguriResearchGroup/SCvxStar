function v_norm = unit(v)
% UNIT - returns the unit vector of v
%   v_norm = unit(v)
%   v - vector
%   v_norm - unit vector of v
%

n = norm(v);
if n == 0
    v_norm = v;
else
    v_norm = v/n;
end

end