%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Adam Hallberg, Oscar Jemsson
% Date:    2025-11-26
% Status:  Incomplete
%
% Comments:
%   Cleaned version of the code for homework assignment 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
close all;
clear;
rng('default');
clc;

% Toggle if you would like to plot or not.
display = true;

% Setup Variables 
A       = [0.5, 1; 0, 0.5];
B       = [0;1];
R       = 1;
S       = eye(2);
gamma   = 0.9;
N       = 10;
epsilon = 1e-14;
L_0 = [0,0];

% phi is taken from example 11.3
phi = @(x,u) [x(1)^2, x(2)^2, u^2, 2*x(1)*x(2), 2*x(1)*u, 2*x(2)*u];


% Set Sample Size. NOTE, has to be \geq 6.
r = 6; 

% The X_history are Just Needed for Plotting.
L_history       = [];
error_history   = [];


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
        break;
    end
    
    L = L_new;
    L_history = [L_history; L];
    k = k+1;
end

if(display)
    plot_convergence(k, L_history, error_history, epsilon);
end
    

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


function plot_convergence(k, L_history, error_history, epsilon)

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
