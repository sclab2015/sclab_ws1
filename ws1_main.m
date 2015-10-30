% main file for worksheet 1 of the scientific computing lab

% definition of analytic solution and initial conditions
analytic_sol = @(t) 10./(1+9*exp(-t));
y0 = 1;
t0 = 0;
n_stepsizes = 4;
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
euler_vals = cell(1,n_stepsizes);
heun_vals = cell(1,n_stepsizes);
rk4_vals = cell(1,n_stepsizes);
% allocate cells to store error values
euler_err_red = zeros(1,n_stepsizes);
heun_err_red = zeros(1,n_stepsizes);
rk4_err_red = zeros(1,n_stepsizes);
euler_err_app = zeros(1,n_stepsizes);
heun_err_app = zeros(1,n_stepsizes);
rk4_err_app = zeros(1,n_stepsizes);

compute_err = @(v_app, v_true, step_size) sqrt(step_size/5*sum((v_app-v_true).^2));

% compute values
analytic_vals = analytic_sol(times);
    
for s = 1:n_stepsizes
    step_size = 2.0^-(s-1);
    rk4_vals{s} = rungekutta4(y0, t0, 5, step_size, @test_func);
    euler_vals{s} = Euler(y0, t0, 5, step_size, @test_func);
    heun_vals{s} = Heun(y0, t0, 5, step_size, @test_func);
end


% compute error tables
for s = 1:n_stepsizes
   step_size = 2.0^-(s-1);
   stride = ceil(8*step_size);
   euler_err_red(s) = sqrt(step_size/5*sum((euler_vals{s}-analytic_vals(1:stride:end)).^2));
   euler_err_app(s) = sqrt(step_size/5*sum((euler_vals{s}-euler_vals{n_stepsizes}(1:stride:end)).^2));
   
   heun_err_red(s) = sqrt(step_size/5*sum((heun_vals{s}-analytic_vals(1:stride:end)).^2));
   heun_err_app(s) = sqrt(step_size/5*sum((heun_vals{s}-heun_vals{n_stepsizes}(1:stride:end)).^2));
   
   rk4_err_red(s) = sqrt(step_size/5*sum((rk4_vals{s}-analytic_vals(1:stride:end)).^2));
   rk4_err_app(s) = sqrt(step_size/5*sum((rk4_vals{s}-rk4_vals{n_stepsizes}(1:stride:end)).^2));
end

% display error values
fprintf('Error values of the Euler method:\n')
for s = 1:n_stepsizes
   fprintf('for %s absolute error is: %1.7g\n', step_strings{s}, euler_err_red(s))
end
for s = 1:n_stepsizes
   fprintf('for %s approximated error is: %1.7g\n', step_strings{s}, euler_err_app(s))
end
fprintf('Error values of the Heun method:\n')
for s = 1:n_stepsizes
   fprintf('for %s absolute error is: %1.7g\n', step_strings{s}, heun_err_red(s))
end
for s = 1:n_stepsizes
   fprintf('for %s approximated error is: %1.7g\n', step_strings{s}, heun_err_app(s))
end
fprintf('Error values of the Runge-Kutta method of fourth order:\n')
for s = 1:n_stepsizes
   fprintf('for %s absolute error is: %1.7g\n', step_strings{s}, rk4_err_red(s))
end
for s = 1:n_stepsizes
   fprintf('for %s approximated error is: %1.7g\n', step_strings{s}, rk4_err_app(s))
end



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
for s = 1:n_stepsizes
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
for s = 1:n_stepsizes
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
for s = 1:n_stepsizes
    rk4_plots(s) = plot(times(1:8*2^-(s-1):end), rk4_vals{s}, colors{s});
end
legend(rk4_plots, step_strings, 'Location', 'NorthWest');
xlabel('time');
ylabel('population');

shg