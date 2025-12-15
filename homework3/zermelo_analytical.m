% Zermelo Problem - Analytical Solution
% Exercise 1.1(d)

clear all;
close all;
clc;

% Parameters
tf = 1;      % Final time
w = 1;       % Ship speed relative to water
X0 = [0; 0]; % Initial position [x; y]

% Time span
tspan = [0 tf];

% Solve ODE with optimal control
[t, X] = ode45(@(t, X) zermelo_dynamics(t, X, w, tf), tspan, X0);

% Extract x and y trajectories
x = X(:, 1);
y = X(:, 2);

% Calculate optimal control trajectory
theta = arrayfun(@(t) atan(tf - t), t);

% Calculate final cost (objective value)
J_opt = -x(end);
fprintf('Numerical optimal objective value: J = %.6f\n', J_opt);
fprintf('Final position: x(tf) = %.6f, y(tf) = %.6f\n', x(end), y(end));

% Plot results
figure('Position', [100 100 1200 400]);

% Subplot 1: Trajectory in x-y plane
subplot(1, 2, 1);
plot(x, y, 'b-', 'LineWidth', 2);
hold on;
plot(x(1), y(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(x(end), y(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on;
xlabel('x');
ylabel('y');
title('Ship Trajectory');
legend('Trajectory', 'Start', 'End', 'Location', 'best');
axis equal;

% Subplot 3: Optimal control vs time
subplot(1, 2, 2);
plot(t, theta * 180/pi, 'k-', 'LineWidth', 2);
grid on;
xlabel('Time t');
ylabel('\theta(t) [degrees]');
title('Optimal Control (Heading Angle)');

% Add visualization of water current
figure('Position', [100 600 800 600]);
% Create a grid for vector field
[X_grid, Y_grid] = meshgrid(linspace(min(x)-0.2, max(x)+0.2, 15), ...
                             linspace(min(y)-0.2, max(y)+0.2, 15));
% Water current velocity: v(y) = y
U_current = Y_grid;  % x-component of water velocity
V_current = zeros(size(Y_grid));  % y-component is zero

% Plot water current as vector field
quiver(X_grid, Y_grid, U_current, V_current, 0.5, 'Color', [0.7 0.7 1], 'LineWidth', 1);
hold on;

% Plot ship trajectory
plot(x, y, 'b-', 'LineWidth', 2.5);
plot(x(1), y(1), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
plot(x(end), y(end), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

% Plot ship heading at several points
n_arrows = 10;
indices = round(linspace(1, length(t), n_arrows));
for i = indices
    % Ship velocity relative to water
    dx_ship = w * cos(theta(i));
    dy_ship = w * sin(theta(i));
    % Scale for visualization
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