% Zermelo Problem - Discretization Method
% Exercise 1.2 - Using fmincon with constrained formulation

clear all;
close all;
clc;

% Parameters
tf = 1;      % Final time
w = 1;       % Ship speed relative to water
Xi = [0; 0]; % Initial position [x; y]

% Discretization parameters
N = 50;      % Number of time steps
T = tf / N;  % Sample time

% Y = [x[0], y[0], theta[0], x[1], y[1], theta[1], ..., x[N], y[N], theta[N]]
% Size: 3*(N+1)
Y0 = zeros(3*(N+1), 1);

% Initial guess using analytical solution
t_guess = linspace(0, tf, N+1)';
for k = 0:N
    Y0(3*k+1) = t_guess(k+1);                    % x[k]
    Y0(3*k+2) = 0.4 * t_guess(k+1);              % y[k]
    Y0(3*k+3) = atan(tf - t_guess(k+1));         % theta[k]
end

% Set up fmincon options
options = optimset('fmincon');
options = optimset(options, 'Algorithm', 'interior-point');
options = optimset(options, 'GradObj', 'on');
options = optimset(options, 'GradConstr', 'on');
options = optimset(options, 'MaxFunEvals', 15000);
options = optimset(options, 'Display', 'iter');
options = optimset(options, 'MaxIter', 1000);

% No linear constraints
Aeq = [];
Beq = [];

fprintf('\nSolving with GRADIENTS enabled...\n');
tic;
[Y_opt, fval] = fmincon(@(Y) zermeloCostFun(Y, N), Y0, [], [], ...
    Aeq, Beq, [], [], @(Y) zermeloNonlcon(Y, Xi, T, w, N), options);
time_with_grad = toc;

fprintf('\n=================================================\n');
fprintf('RESULTS WITH GRADIENTS:\n');
fprintf('=================================================\n');
fprintf('  Optimal cost:           J = %.6f\n', fval);
fprintf('  Final position:         x(tf) = %.6f\n', Y_opt(3*N+1));
fprintf('                          y(tf) = %.6f\n', Y_opt(3*N+2));
fprintf('  Computation time:       %.3f seconds\n', time_with_grad);
fprintf('=================================================\n');

% Extract solution
x_opt = Y_opt(1:3:end);
y_opt = Y_opt(2:3:end);
theta_opt = Y_opt(3:3:end);
t_opt = linspace(0, tf, N+1)';

% Verify constraints
[c, ceq] = zermeloNonlcon(Y_opt, Xi, T, w, N);
fprintf('\nConstraint Verification:\n');
fprintf('  Max constraint violation: %.6e\n', max(abs(ceq)));

% Compare with analytical solution
theta_analytical = atan(tf - t_opt);

fprintf('\nComparison with Analytical Solution:\n');
fprintf('  Analytical solution:    J = -1.147794\n');
fprintf('  Error:                  %.6e (%.2f%%)\n', ...
    abs(fval + 1.147794), 100*abs(fval + 1.147794)/1.147794);

% Now solve WITHOUT gradients for comparison (Exercise 1.2d)
fprintf('\n=================================================\n');
fprintf('Solving WITHOUT gradients...\n');
fprintf('=================================================\n');

options_no_grad = optimset(options, 'GradObj', 'off', 'GradConstr', 'off');
%options_no_grad = optimset(options_no_grad, 'Display', 'iter');

tic;
[Y_opt_no_grad, fval_no_grad] = fmincon(@(Y) zermeloCostFun(Y, N), Y0, [], [], ...
    Aeq, Beq, [], [], @(Y) zermeloNonlcon(Y, Xi, T, w, N), options_no_grad);
time_no_grad = toc;

fprintf('\n=================================================\n');
fprintf('RESULTS WITHOUT GRADIENTS:\n');
fprintf('=================================================\n');
fprintf('  Optimal cost:           J = %.6f\n', fval_no_grad);
fprintf('  Computation time:       %.3f seconds\n', time_no_grad);
fprintf('=================================================\n');

fprintf('\n=================================================\n');
fprintf('COMPARISON (Exercise 1.2d):\n');
fprintf('=================================================\n');
fprintf('Time:\n');
fprintf('  With gradients:         %.3f seconds\n', time_with_grad);
fprintf('  Without gradients:      %.3f seconds\n', time_no_grad);
if time_with_grad < time_no_grad
    fprintf('  Result:                 %.2fx FASTER with gradients\n', time_no_grad / time_with_grad);
else
    fprintf('  Result:                 %.2fx SLOWER with gradients\n', time_with_grad / time_no_grad);
end
fprintf('\nAccuracy:\n');
fprintf('  Analytical solution:    J = -1.147794\n');
fprintf('  Error (with grad):      %.6e\n', abs(fval + 1.147794));
fprintf('  Error (no grad):        %.6e\n', abs(fval_no_grad + 1.147794));
fprintf('  Cost difference:        %.6e\n', abs(fval - fval_no_grad));
fprintf('\nConclusion:\n');
if time_with_grad < time_no_grad
    fprintf('  Analytical gradients are BENEFICIAL:\n');
    fprintf('  - Faster (%.2fx speedup)\n', time_no_grad / time_with_grad);
    fprintf('  - Similar or better accuracy\n');
else
    fprintf('  For this SMALL problem (N=%d), gradients are SLOWER:\n', N);
    fprintf('  - Gradient computation overhead (%.2fx slower)\n', time_with_grad / time_no_grad);
    fprintf('  - Accuracy is similar either way\n');
    fprintf('  - For LARGER problems, gradients typically help!\n');
    fprintf('  - Finite differences may be adequate for small N\n');
end
fprintf('=================================================\n');

% Create plots
if ~exist('figures', 'dir')
    mkdir('figures');
end

figure('Position', [100 100 1400 900]);

% Subplot 1: Trajectory
subplot(2, 3, 1);
plot(x_opt, y_opt, 'b-', 'LineWidth', 2);
hold on;
plot(x_opt(1), y_opt(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(x_opt(end), y_opt(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on;
xlabel('x');
ylabel('y');
title('Optimal Trajectory');
legend('Trajectory', 'Start', 'End', 'Location', 'best');
axis equal;

% Subplot 2: States vs time
subplot(2, 3, 2);
plot(t_opt, x_opt, 'b-', 'LineWidth', 2);
hold on;
plot(t_opt, y_opt, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time t');
ylabel('Position');
title('States vs Time');
legend('x(t)', 'y(t)', 'Location', 'best');

% Subplot 3: Control comparison
subplot(2, 3, 3);
plot(t_opt, theta_opt * 180/pi, 'b-', 'LineWidth', 2);
hold on;
plot(t_opt, theta_analytical * 180/pi, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time t');
ylabel('\theta(t) [degrees]');
title('Optimal Control');
legend('Discretization', 'Analytical', 'Location', 'best');

% Subplot 4: Control error
subplot(2, 3, 4);
plot(t_opt, (theta_opt - theta_analytical)*180/pi, 'k-', 'LineWidth', 2);
grid on;
xlabel('Time t');
ylabel('Error [degrees]');
title('Control Error');

% Subplot 5: Position error
subplot(2, 3, 5);
[~, X_analytical] = ode45(@(t, X) zermelo_dynamics(t, X, w, tf), t_opt, Xi);
pos_error = sqrt((x_opt - X_analytical(:,1)).^2 + (y_opt - X_analytical(:,2)).^2);
plot(t_opt, pos_error, 'k-', 'LineWidth', 2);
grid on;
xlabel('Time t');
ylabel('Position Error');
title('Position Error');

% Subplot 6: Constraint satisfaction
subplot(2, 3, 6);
dynamics_error = zeros(N, 1);
for k = 1:N
    x_pred = x_opt(k) + T * (w * cos(theta_opt(k)) + y_opt(k));
    y_pred = y_opt(k) + T * w * sin(theta_opt(k));
    dynamics_error(k) = sqrt((x_opt(k+1) - x_pred)^2 + (y_opt(k+1) - y_pred)^2);
end
semilogy(t_opt(1:N), dynamics_error, 'k-', 'LineWidth', 2);
grid on;
xlabel('Time t');
ylabel('Constraint Error');
title('Dynamics Constraint Satisfaction');

print('-dpng', '-r300', 'figures/zermelo_fmincon_results.png');
fprintf('\nSaved: figures/zermelo_fmincon_results.png\n');



%% Helper function: System dynamics
function dXdt = zermelo_dynamics(t, X, w, tf)
    % X = [x; y]
    % Optimal control law
    theta_opt = atan(tf - t);
    
    % Water current velocity (linear with y)
    v_y = X(2);
    
    % State derivatives
    dx_dt = w * cos(theta_opt) + v_y;
    dy_dt = w * sin(theta_opt);
    
    dXdt = [dx_dt; dy_dt];
end