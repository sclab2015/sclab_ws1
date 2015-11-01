function [ y ] = Euler(y0, t0, t_end, dt, f)
% TODO: Documentation.

N = ceil((t_end - t0) / dt);

y = zeros(1, N + 1);
y(1) = y0;

% TODO: Wenn t_end - t0 nicht durch dt teilbar ist gehen wir am Ende zu weit!
for i = 1:N
    t = t0 + (i - 1) * dt;
    k1 = dt * f(t, y(i));
    y(i + 1) = y(i) + k1;
end

end
