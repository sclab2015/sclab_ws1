% main file for worksheet 1 of the scientific computing lab

% definition of the right-hand side, initial conditions and the analytic
% solution
f = @(t, y) (1 - t / 10) * t;
analytic_sol = @(t) 10 ./ (1 + 9 * exp(-t));
y0 = 1;
t0 = 0;
t_end = 5;

% definition of time steps and colors for plotting
steps  = 2.^(0:-1:-3);  % [1; .5; .25; .125]
colors = {'r+', 'g+', 'c+', 'm+'};
assert(length(steps) == length(colors), 'number of steps ~= number of colors');
times = t0 : steps(end) : t_end;

% allocate memory for plot handles and define legend entries
euler_plots = gobjects(1, 3);
heun_plots = gobjects(1, 3);
rk4_plots = gobjects(1, 3);
step_strings = {'step 1', 'step 1/2', 'step 1/4', 'step 1/8'};
assert(length(steps) == length(step_strings), ...
       'number of steps ~= number of labels');

% allocate cells to store computed function values for different step sizes
euler_vals = cell(1, length(steps));
heun_vals = cell(1, length(steps));
rk4_vals = cell(1, length(steps));

% allocate cells to store error values
euler_err = zeros(1, length(steps));  % TODO: _red ist wohl falsch!
heun_err = zeros(1, length(steps));
rk4_err = zeros(1, length(steps));
euler_err_app = zeros(1, length(steps));
heun_err_app = zeros(1, length(steps));
rk4_err_app = zeros(1, length(steps));

compute_err = @(p, p_better, dt) sqrt(dt / 5 * sum((p - p_better).^2));

% compute analytic values
analytic_vals = analytic_sol(times);

% compute numerical approximation
for i = 1:length(steps)
    dt = steps(i);
    euler_vals{i} = Euler(y0, t0, t_end, dt, f);
    heun_vals{i} = Heun(y0, t0, t_end, dt, f);
    rk4_vals{i} = RK4(y0, t0, t_end, dt, f);
end

% compute error tables
for i = 1:length(steps)
    dt = steps(i);
    % based on the assumption that we use stepsizes 2^-k
    stride = ceil(dt / steps(end));

    euler_err(i) = compute_err(euler_vals{i}, ...
                               analytic_vals(1:stride:end), dt);
    euler_err_app(i) = compute_err(euler_vals{i}, ...
                                   euler_vals{length(steps)}(1:stride:end), dt);

    heun_err(i) = compute_err(heun_vals{i}, analytic_vals(1:stride:end), ...
                              dt);
    heun_err_app(i) = compute_err(heun_vals{i}, ...
                                  heun_vals{length(steps)}(1:stride:end), dt);

    rk4_err(i) = compute_err(rk4_vals{i}, analytic_vals(1:stride:end), dt);
    rk4_err_app(i) = compute_err(rk4_vals{i}, ...
                                 rk4_vals{length(steps)}(1:stride:end), dt);
end

% TODO: Implement error reduction here!

% print error table
% euler
fprintf('explicit Euler method (q = 1):\n');
fprintf('dt\t\t');
for dt = steps
    fprintf('%f\t', dt);
end
fprintf('\n');
fprintf('error\t\t');
for e = euler_err
    fprintf('%f\t', e);
end
fprintf('\n')
%fprintf('error red.\t')
% TODO: Implement error reduction.
%for e = euler_err_red
%    fprintf('%f\t', e);
%end
%fprintf('\n')
fprintf('error app.\t')
for e = euler_err_app
    fprintf('%f\t', e);
end
fprintf('\n\n')

% heun
fprintf('method of Heun (q = 2):\n');
fprintf('dt\t\t');
for dt = steps
    fprintf('%f\t', dt);
end
fprintf('\n');
fprintf('error\t\t');
for e = heun_err
    fprintf('%f\t', e);
end
fprintf('\n')
%fprintf('error red.\t')
% TODO: Implement error reduction.
%for e = euler_err_red
%    fprintf('%f\t', e);
%end
%fprintf('\n')
fprintf('error app.\t')
for e = heun_err_app
    fprintf('%f\t', e);
end
fprintf('\n\n')

% euler
fprintf('Runge-Kutta method (q = 4):\n');
fprintf('dt\t\t');
for dt = steps
    fprintf('%f\t', dt);
end
fprintf('\n');
fprintf('error\t\t');
for e = rk4_err
    fprintf('%f\t', e);
end
fprintf('\n')
%fprintf('error red.\t')
% TODO: Implement error reduction.
%for e = euler_err_red
%    fprintf('%f\t', e);
%end
%fprintf('\n')
fprintf('error app.\t')
for e = rk4_err_app
    fprintf('%f\t', e);
end
fprintf('\n\n')

% create figure (fullscreen) with a subplot for each method
figure('units','normalized','outerposition',[0 0 1 1])
hold on

% generate plot for euler method
subplot(1, 3, 1)
hold on
title('Integration With Euler''s Method')
% plot analytic solution
plot(times, analytic_vals, 'b');
% plot different timesteps
for i = 1:length(steps)
    euler_plots(i) = plot(times(1:8*2^ - (i - 1):end), euler_vals{i}, ...
                          colors{i});
end
legend(euler_plots, step_strings, 'Location', 'NorthWest');
xlabel('time');
ylabel('population');

% generate plot for heun method
subplot(1, 3, 2)
hold on
title('Integration With Heun''s Method')
% plot analytic solution
plot(times, analytic_vals, 'b');
% plot different timesteps
for i = 1:length(steps)
    heun_plots(i) = plot(times(1:8*2^ - (i - 1):end), heun_vals{i}, colors{i});
end
legend(heun_plots, step_strings, 'Location', 'NorthWest');
xlabel('time');
ylabel('population');

% generate plot for 4th order runge kutta
subplot(1, 3, 3)
hold on
title('Integration With 4th Order Runge-Kutta')
% plot analytic solution
plot(times, analytic_vals, 'b');
% plot different timesteps
for i = 1:length(steps)
    rk4_plots(i) = plot(times(1:8*2^ - (i - 1):end), rk4_vals{i}, colors{i});
end
legend(rk4_plots, step_strings, 'Location', 'NorthWest');
xlabel('time');
ylabel('population');

shg