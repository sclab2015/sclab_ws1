% main file for worksheet 1 of the scientific computing lab

% definition of analytic solution and initial conditions
analytic_sol = @(t) 10./(1+9*exp(-t));
y0 = 1;
t0 = 0;
% definition of time steps and colors for plotting
steps  = [1.0, 0.5, 0.25, 0.125];
colors = {'r+', 'g+', 'c+', 'm+'};
times = linspace(0, 5, 41);
% allocate memory for plot handles and define legend entries
euler_plots = gobjects(1,3);
heun_plots = gobjects(1,3);
rk4_plots = gobjects(1,3);
step_strings = {'step 1', 'step 1/2', 'step 1/4', 'step1/8'};
% allocate cells to store computed function values for different step sizes
euler_vals = cell(1,4);
heun_vals = cell(1,4);
rk4_vals = cell(1,4);


% compute values
analytic_vals = analytic_sol(times);
    
for s = 1:4
    step_size = 2.0^-(s-1);
    rk4_vals{s} = rungekutta4(y0, t0, 5, step_size, @test_func);
    euler_vals{s} = Euler(y0, t0, 5, step_size, @test_func);
    heun_vals{s} = Heun(y0, t0, 5, step_size, @test_func);
end


% compute error tables





% create figure with a subplot for each method
figure()
hold on

% generate plot for euler method
subplot(1,3,1)
hold on
title('Integration With Euler''s Method')
% plot analytic solution
plot(times, analytic_vals, 'b');
% plot different timesteps
for s = 1:4
    euler_plots(s) = plot(times(1:8*2^-(s-1):end), euler_vals{s}, colors{s});
end
legend(euler_plots, step_strings, 'Location', 'NorthWest');
xlabel('time');
ylabel('population');

% generate plot for heun method
subplot(1,3,2)
hold on
title('Integration With Heun''s Method')
% plot analytic solution
plot(times, analytic_vals, 'b');
% plot different timesteps
for s = 1:4
    heun_plots(s) = plot(times(1:8*2^-(s-1):end), heun_vals{s}, colors{s});
end
legend(heun_plots, step_strings, 'Location', 'NorthWest');
xlabel('time');
ylabel('population');

% generate plot for 4th order runge kutta
subplot(1,3,3)
hold on
title('Integration With 4th Order Runge-Kutta')
% plot analytic solution
plot(times, analytic_vals, 'b');
% plot different timesteps
for s = 1:4
    rk4_plots(s) = plot(times(1:8*2^-(s-1):end), rk4_vals{s}, colors{s});
end
legend(rk4_plots, step_strings, 'Location', 'NorthWest');
xlabel('time');
ylabel('population');

shg