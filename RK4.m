function [ y ] = RK4(y0, t0, t_end, dt, f)
  % TODO: Documentation.

N = ceil((t_end - t0) / dt);

y = zeros(1, N + 1);
y(1) = y0;

% TODO: Wenn t_end - t0 nicht durch dt teilbar ist gehen wir am Ende zu weit!
for i = 1:N
    t = t0 + (i - 1) * dt;
    k1 = dt * f(t, y(i));
    k2 = dt * f(t + .5 * dt, y(i) + .5 * k1);
    k3 = dt * f(t + .5 * dt, y(i) + .5 * k2);
    k4 = dt * f(t + dt, y(i) + k3);
    y(i + 1) = y(i) + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4);
end

end