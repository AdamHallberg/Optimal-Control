% Zermelo Problem - Gradient Method
% Exercise 1.3 - Using gradient search approach

clear all;
close all;
clc;

tf = 1;      % Final time
w = 1;       % Ship speed relative to water
Xi = [0; 0];

max_iter = 1000;
tol = 1e-6;
alpha = 0.8;

N_time = 100;
t_span = linspace(0, tf, N_time);

theta0 = zeros(N_time, 1);

fprintf('\nStarting gradient descent...\n');
fprintf('Initial step size alpha = %.3f\n', alpha);

theta = theta0;
J_history = zeros(max_iter, 1);
gradient_norms = zeros(max_iter, 1);

tic;
for iter = 1:max_iter
    % solve state equations
    [t, X] = ode45(@(t, X) zermelo_state_dynamics(t, X, theta, t_span, w), ...
                   t_span, Xi);
    x = X(:, 1);
    y = X(:, 2);
    
    % Compute cost
    J = -x(end);
    J_history(iter) = J;
    
    % solve adjoint equations (Backwards)
    lambda_f = [-1; 0];
    
    [t_back, Lambda] = ode45(@(t, lambda) zermelo_adjoint_dynamics(t, lambda, x, y, theta, t_span, w), ...
                             flip(t_span), lambda_f);
    
    % Reverse to get forward time
    Lambda = flip(Lambda, 1);
    lambda_x = Lambda(:, 1);
    lambda_y = Lambda(:, 2);
    
    % dH/dtheta = -lambda_x * w * sin(theta) + lambda_y * w * cos(theta)
    grad_theta = zeros(N_time, 1);
    for i = 1:N_time
        grad_theta(i) = -lambda_x(i) * w * sin(theta(i)) + ...
                        lambda_y(i) * w * cos(theta(i));
    end
    
    % Check convergence
    grad_norm = norm(grad_theta);
    gradient_norms(iter) = grad_norm;
    
    if mod(iter, 50) == 0 || iter == 1
        fprintf('Iter %4d: J = %.6f, ||grad|| = %.6e\n', iter, J, grad_norm);
    end
    
    if grad_norm < tol
        fprintf('\nConverged at iteration %d\n', iter);
        break;
    end
    
    % Gradient descent update
    theta = theta - alpha * grad_theta;
end
time_elapsed = toc;

[t_final, X_final] = ode45(@(t, X) zermelo_state_dynamics(t, X, theta, t_span, w), ...
                           t_span, Xi);
x_final = X_final(:, 1);
y_final = X_final(:, 2);
J_final = -x_final(end);

fprintf('\n=================================================\n');
fprintf('RESULTS:\n');
fprintf('=================================================\n');
fprintf('  Optimal cost:           J = %.6f\n', J_final);
fprintf('  Final position:         x(tf) = %.6f\n', x_final(end));
fprintf('                          y(tf) = %.6f\n', y_final(end));
fprintf('  Iterations:             %d\n', iter);
fprintf('  Computation time:       %.3f seconds\n', time_elapsed);
fprintf('  Final gradient norm:    %.6e\n', grad_norm);
fprintf('=================================================\n');

% Compare with analytical and discretization solutions
theta_analytical = atan(tf - t_span');

fprintf('\nComparison:\n');
fprintf('  Analytical solution:    J = -1.147794\n');
fprintf('  Gradient method:        J = %.6f\n', J_final);
fprintf('  Error:                  %.6e (%.2f%%)\n', ...
    abs(J_final + 1.147794), 100*abs(J_final + 1.147794)/1.147794);
fprintf('=================================================\n');

if ~exist('figures', 'dir')
    mkdir('figures');
end

figure('Position', [100 100 1400 900]);

subplot(2, 2, 1);
plot(x_final, y_final, 'b-', 'LineWidth', 2);
hold on;
plot(x_final(1), y_final(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(x_final(end), y_final(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on;
xlabel('x');
ylabel('y');
title('Optimal Trajectory (Gradient Method)');
legend('Trajectory', 'Start', 'End', 'Location', 'best');
axis equal;

subplot(2, 2, 2);
plot(t_span, theta * 180/pi, 'b-', 'LineWidth', 2);
hold on;
plot(t_span, theta_analytical * 180/pi, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time t');
ylabel('\theta(t) [degrees]');
title('Optimal Control');
legend('Gradient Method', 'Analytical', 'Location', 'best');

subplot(2, 2, 3);
plot(1:iter, J_history(1:iter), 'k-', 'LineWidth', 2);
grid on;
xlabel('Iteration');
ylabel('Cost J');
title('Convergence History');

subplot(2, 2, 4);
semilogy(1:iter, gradient_norms(1:iter), 'k-', 'LineWidth', 2);
grid on;
xlabel('Iteration');
ylabel('||∇J||');
title('Gradient Norm History');

print('-dpng', '-r300', 'figures/zermelo_gradient_method.png');
fprintf('\nSaved: figures/zermelo_gradient_method.png\n');

