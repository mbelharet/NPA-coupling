function d_len = get_day_length(Time,Lat)
p = 0.833;
theta = 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (Time + 1 - 186))); % eq. 1 % (length(t)xlength(L))
phi = asin(0.39795 .* cos(theta));    % eq. 2  % (length(t)xlength(L))
a = (sin(p .* pi ./ 180.) + sin(Lat .* pi ./ 180.) .* sin(phi)) ./ (cos(Lat * pi / 180.) .* cos(phi));
a = min(1, a);
a = max(-1, a);
d_len = 1 - (1 ./ pi) .* acos(a);
% d_len = min(1, d_len);   % prevents values > 1
% d_len = max(1/24, d_len);

end