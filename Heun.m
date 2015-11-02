function [ y ] = Heun(y0, t0, t_end, dt, f)
%Heun Numerically solves a scalar differential equation of the form
%     y'(t) = f(t, y(t)), y(t0) = y0 using Heun's method (second-order) with
%     constant stepsize dt on the interval [t0, t_end].

if mod(t_end - t0, dt)
    fprintf('Warning: Integration is done until a full multiple of time step dt is reached!');
end

N = ceil((t_end - t0) / dt);

y = zeros(1, N + 1);
y(1) = y0;

% TODO: Wenn t_end - t0 nicht durch dt teilbar ist gehen wir am Ende zu weit!
for i = 1:N
    t = t0 + (i - 1) * dt;
    k1 = dt * f(t, y(i));
    k2 = dt * f(t + dt, y(i) + k1);
    y(i + 1) = y(i) + 1/2 * (k1 + k2);
end

end
