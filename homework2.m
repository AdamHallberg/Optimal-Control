clear
clc
A = [0.5, 1; 0, 0.5];
B = [0;1];
R = 1;
S = eye(2);
gamma = 0.9;
N = 10;
epsilon = 1e-14;
phi = @(x,u) [x(1)^2, x(2)^2, u^2, 2*x(1)*x(2), 2*x(1)*u, 2*x(2)*u];
L_0 = [0,0];
r = 6; % At least 6, (sample size)

% Storage for visualization
L_history = [];
error_history = [];

L = L_0;
L_history = [L_history; L];
k = 1;

% Create figure for real-time visualization
figure('Position', [100, 100, 1200, 800]);

while (true)
    x = randn(2,r);
    u = randn(1,r);
    beta = zeros(r,1);
    phi_mat = zeros(r,6);
    
    for s = 1:r
        beta(s) = get_beta(x(:,s), u(:,s), N, L, S, R, gamma, A, B);
        phi_mat(s,:) = phi(x(:,s), u(:,s));
    end
    
    % Phi'Phi a = Phi' Beta
    a = (phi_mat'*phi_mat)\(phi_mat'*beta); % LS solution
    P = [a(1), a(4); a(4), a(2)];
    r_tilde = [a(5); a(6)];
    q_tilde = a(3);
    L_new = q_tilde\r_tilde';
    
    policy_change = norm(L - L_new);
    error_history = [error_history; policy_change];
    
    if (policy_change < epsilon)
        L_history = [L_history; L_new];
        break;
    end
    
    L = L_new;
    L_history = [L_history; L];
    k = k+1;
end
% Visualization
% Subplot 1: Policy convergence (L gains)
subplot(1, 2, 1);
plot(0:k, L_history(:,1), 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(0:k, L_history(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('Iteration', 'FontSize', 12);
ylabel('Gain Value', 'FontSize', 12);
title('Policy Convergence', 'FontSize', 14, 'FontWeight', 'bold');
legend('L(1)', 'L(2)', 'Location', 'best');

% Subplot 2: Policy change (error)
subplot(1, 2, 2);
semilogy(1:length(error_history), error_history, 'm-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
yline(epsilon, 'k--', 'LineWidth', 1.5, 'Label', '\epsilon threshold');
grid on;
xlabel('Iteration', 'FontSize', 12);
ylabel('||L_{k+1} - L_k||', 'FontSize', 12);
title('Policy Change per Iteration', 'FontSize', 14, 'FontWeight', 'bold');


function beta = get_beta(x,u,N,L, S, R,gamma,A,B) 
%  x, u  column vectors
    beta = x'*S*x + u'*R*u;
    x = A*x + B*u;
    for i = 1:N-1
        beta = beta + gamma^i * (x'*S*x + (-L*x)'*R*(-L*x));
        x = (A - B*L)*x;
    end
end