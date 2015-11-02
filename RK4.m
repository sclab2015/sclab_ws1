function [ y ] = RK4(y0, t0, t_end, dt, f)
%RK4 Numerically solves a scalar differential equation of the form
%    y'(t) = f(t, y(t)), y(t0) = y0 using the fourth-order Runge-Kutta scheme
%    with constant stepsize dt on the interval [t0, t_end].

if mod(t_end - t0, dt)
    fprintf('Warning: Integration is done until a full multiple of time step dt is reached!');
end

N = ceil((t_end - t0) / dt); %lalala

y = zeros(1, N + 1);
y(1) = y0;

for i = 1:N
    t = t0 + (i - 1) * dt;
    k1 = dt * f(t, y(i));
    k2 = dt * f(t + .5 * dt, y(i) + .5 * k1);
    k3 = dt * f(t + .5 * dt, y(i) + .5 * k2);
    k4 = dt * f(t + dt, y(i) + k3);
    y(i + 1) = y(i) + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4);
end

end
