%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Adam Hallberg, Oscar Jemsson
% Date:    2025-11-26
% Status:  Incomplete
%
% Comments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
close all;
clear;
%rng('default');
clc;

display = true;

% Setup variables
N_base = 15;
r_base = 10; 
simulations = 40;
epsilon = 1e-10;


r_values      = [];
r_mean_values = [];
for r = 6:15
% Run simulation
[iterations, L_history, error_history] = simulation(r, N_base, simulations, epsilon);

r_values = [r_values; r];
r_mean_values = [r_mean_values; mean(iterations);];
end

N_values      = [];
N_mean_values = [];
for N = 3:30
% Run simulation
[iterations, L_history, error_history] = simulation(r_base, N, simulations, epsilon);

N_values = [N_values; N];
N_mean_values = [N_mean_values; mean(iterations);];
end


if(display)
    % Display how r affects mean convergence
    r_plot(r_values, r_mean_values)

    % Dsiplay How N affects mean convergence
    N_plot(N_values, N_mean_values)

end

%%




%% ===================FUNCTIONS============================================
function beta = get_beta(x,u,N,L, S, R,gamma,A,B) 
%  x, u  column vectors
    beta = x'*S*x + u'*R*u;
    x = A*x + B*u;
    for i = 1:N-1
        beta = beta + gamma^i * (x'*S*x + (-L*x)'*R*(-L*x));
        x = (A - B*L)*x;
    end
end

function r_plot(r_values, r_mean_values)
    figure('Position',[200 200 900 600]);
    
    % Plot mean iterations
    plot(r_values, r_mean_values, 'o-', ...
        'LineWidth', 1.8, 'MarkerSize', 6);
    hold on;
    
    
    grid on;
    
    xlabel('Number of samples r','FontSize',14);
    ylabel('Mean iterations to convergence','FontSize',14);
    title('Effect of Sample Size on Convergence Speed','FontSize',16,'FontWeight','bold');
    
    legend('Mean iterations','Location','northeast');
    
    set(gca,'FontSize',12);

end

function N_plot(N_values, N_mean_values)
    figure('Position',[200 200 900 600]);
    
    % Plot mean iterations
    plot(N_values, N_mean_values, 'o-', ...
        'LineWidth', 1.8, 'MarkerSize', 6);
    hold on;
    
    
    grid on;
    
    xlabel('N','FontSize',14);
    ylabel('Mean iterations to convergence','FontSize',14);
    title('Effect of N on Convergence Speed','FontSize',16,'FontWeight','bold');
    
    legend('Mean iterations','Location','northeast');
    
    set(gca,'FontSize',12);

end


function plot_iterations(iterations, r_start)

    % Keep only columns that were actually filled
    valid = any(iterations > 0, 1);
    iterations = iterations(:, valid);

    % Compute statistics
    mean_iters = mean(iterations,1);

    % X-axis: number of samples r
    r_values = r_start:(r_start + length(mean_iters) - 1);

    % Print table
    fprintf("\n===== Convergence vs Sample Size =====\n");
    fprintf("r\tmean iters\n");
    for j = 1:length(r_values)
        fprintf("%d\t%.2f\n", r_values(j), mean_iters(j));
    end

    % Plot
    figure;
    plot(r_values, mean_iters, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
    grid on;
    xlabel('Number of samples r', 'FontSize', 12);
    ylabel('Mean iterations to convergence', 'FontSize', 12);
    title('Effect of Sample Size on Convergence Speed', ...
          'FontSize', 14, 'FontWeight', 'bold');

end




function plot_convergence(L_history, error_history, epsilon)

    k = length(L_history(:,1))-1;
    figure('Position', [100, 100, 1200, 800]);
    
    % Subplot 1: Policy convergence
    subplot(1, 2, 1);
    plot(0:k, L_history(:,1), 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    hold on;
    plot(0:k, L_history(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
    grid on;
    xlabel('Iteration', 'FontSize', 12);
    ylabel('Gain Value', 'FontSize', 12);
    title('Policy Convergence', 'FontSize', 14, 'FontWeight', 'bold');
    legend('L(1)', 'L(2)', 'Location', 'best');
    
    % Subplot 2: Policy change 
    subplot(1, 2, 2);
    semilogy(1:length(error_history), error_history, 'm-o', 'LineWidth', 2, 'MarkerSize', 6);
    hold on;
    yline(epsilon, 'k--', 'LineWidth', 1.5, 'Label', '\epsilon threshold');
    grid on;
    xlabel('Iteration', 'FontSize', 12);
    ylabel('||L_{k+1} - L_k||', 'FontSize', 12);
    title('Policy Change per Iteration', 'FontSize', 14, 'FontWeight', 'bold');
end




function [iterations, L_history, error_history] = simulation(r, N, simulations, epsilon)
    
    % Setup Variables 
    A = [0.5, 1; 0, 0.5];
    B = [0;1];
    R = 1;
    S = eye(2);
    gamma = 0.9;
    phi = @(x,u) [x(1)^2, x(2)^2, u^2, 2*x(1)*x(2), 2*x(1)*u, 2*x(2)*u];
    L_0 = [0,0];
    
        
    iterations = zeros(simulations, 1);
    
    for i = 1:simulations
        % Storage for visualization
        L_history = [];
        error_history = [];
        
        L = L_0;
        L_history = [L_history; L];
        k = 1;
        while (true)
            % Sample x and u
            x = randn(2,r);
            u = randn(1,r);
        
            % Allocate beta and the phi matrix
            beta    = zeros(r,1);
            phi_mat = zeros(r,6); % Taken from 11.3
            
            % This is also taken from 11.3
            for s = 1:r
                beta(s) = get_beta(x(:,s), u(:,s), N, L, S, R, gamma, A, B);
                phi_mat(s,:) = phi(x(:,s), u(:,s));
            end
            
            % We have, Phi'Phi a = Phi' Beta, here ncol < nrow => LS solution
            a = (phi_mat'*phi_mat)\(phi_mat'*beta); 
        
            % With the a matrix we get P, r_tilde and q_tilde according to 11.3
            P = [a(1), a(4); a(4), a(2)];
            r_tilde = [a(5); a(6)];
            q_tilde = a(3);
            L_new = q_tilde\r_tilde';
            
            % For the plot
            policy_change = norm(L - L_new);
            error_history = [error_history; policy_change];
            
            % Terminal condition
            if (policy_change < epsilon)
                L_history = [L_history; L_new];
                iterations(i) = k;
                break;
            end
            
            L = L_new;
            L_history = [L_history; L];
            k = k+1;
        end
    end

    valid = any(iterations > 0, 1);
    iterations = iterations(:, valid);

end


