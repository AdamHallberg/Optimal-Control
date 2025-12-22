% Zermelo Problem - Analytical Solution
% Exercise 1.1(d)

clear all;
close all;
clc;

tf = 1;      % Final time
w = 1;       % Ship speed relative to water
X0 = [0; 0];
tspan = [0 tf];

[t, X] = ode45(@(t, X) zermelo_dynamics(t, X, w, tf), tspan, X0);

x = X(:, 1);
y = X(:, 2);

% optimal control
theta = arrayfun(@(t) atan(tf - t), t);

% objective value
J_opt = -x(end);
fprintf('Numerical optimal objective value: J = %.6f\n', J_opt);
fprintf('Final position: x(tf) = %.6f, y(tf) = %.6f\n', x(end), y(end));


% Add visualization of water current
figure('Position', [100 600 800 600]);
[X_grid, Y_grid] = meshgrid(linspace(min(x)-0.2, max(x)+0.2, 15), ...
                             linspace(min(y)-0.2, max(y)+0.2, 15));
U_current = Y_grid;
V_current = zeros(size(Y_grid));
quiver(X_grid, Y_grid, U_current, V_current, 0.5, 'Color', [0.7 0.7 1], 'LineWidth', 1);
hold on;

% ship trajectory
plot(x, y, 'b-', 'LineWidth', 2.5);
plot(x(1), y(1), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
plot(x(end), y(end), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

n_arrows = 10;
indices = round(linspace(1, length(t), n_arrows));
for i = indices
    dx_ship = w * cos(theta(i));
    dy_ship = w * sin(theta(i));
    scale = 0.15;
    quiver(x(i), y(i), dx_ship*scale, dy_ship*scale, 0, ...
           'Color', 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
end

grid on;
xlabel('x');
ylabel('y');
title('Zermelo Problem: Ship Trajectory with Water Current');
legend('Water current', 'Ship trajectory', 'Start', 'End', 'Ship heading', ...
       'Location', 'best');
axis equal;

fprintf('\nAnalytical solution implemented successfully!\n');

%% System dynamics
function dXdt = zermelo_dynamics(t, X, w, tf)
    theta_opt = atan(tf - t);
    
    v_y = X(2);

    dx_dt = w * cos(theta_opt) + v_y;
    dy_dt = w * sin(theta_opt);
    
    dXdt = [dx_dt; dy_dt];
end