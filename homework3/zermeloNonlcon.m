function [c, ceq, gradc, gradceq] = zermeloNonlcon(Y, Xi, T, w, N)
    % Nonlinear constraints for Zermelo problem
    %
    % Inputs:
    %   Y  - optimization variable [X[0]^T, theta[0], ..., X[N]^T, theta[N]]
    %   Xi - initial state [x0; y0]
    %   T  - sample time
    %   w  - ship speed
    %   N  - number of time steps
    %
    % Outputs:
    %   c       - inequality constraints (none)
    %   ceq     - equality constraints H(Y) = 0
    %   gradc   - gradient of inequality constraints
    %   gradceq - gradient of equality constraints
    
    % No inequality constraints
    c = [];
    
    % Equality constraints: H(Y) = 0
    % H has N+1 constraints (2 for initial condition, then 2 for each step)
    ceq = zeros(2*(N+1), 1);
    
    % Initial condition: X[0] = Xi
    X0 = [Y(1); Y(2)];
    ceq(1:2) = X0 - Xi;
    
    % Dynamic constraints: X[k+1] = F_bar(X[k], theta[k]) for k=0,...,N-1
    for k = 0:N-1
        % Extract X[k] and theta[k]
        x_k = Y(3*k + 1);
        y_k = Y(3*k + 2);
        theta_k = Y(3*k + 3);
        
        % Extract X[k+1]
        x_k1 = Y(3*(k+1) + 1);
        y_k1 = Y(3*(k+1) + 2);
        
        % Discrete-time dynamics using Euler approximation
        % x[k+1] = x[k] + T * (w*cos(theta[k]) + y[k])
        % y[k+1] = y[k] + T * w*sin(theta[k])
        x_k1_pred = x_k + T * (w * cos(theta_k) + y_k);
        y_k1_pred = y_k + T * w * sin(theta_k);
        
        % Constraint: X[k+1] - F_bar(X[k], theta[k]) = 0
        ceq(2*(k+1) + 1) = x_k1 - x_k1_pred;
        ceq(2*(k+1) + 2) = y_k1 - y_k1_pred;
    end
    
    % Compute gradient if requested
    if nargout > 2
        gradc = [];
        
        % gradceq is a matrix of size length(Y) x length(ceq)
        % Each column is the gradient of one constraint
        gradceq = zeros(length(Y), length(ceq));
        
        % Gradient of initial condition constraints
        % ceq(1) = Y(1) - Xi(1), so d/dY(1) = 1
        gradceq(1, 1) = 1;
        % ceq(2) = Y(2) - Xi(2), so d/dY(2) = 1
        gradceq(2, 2) = 1;
        
        % Gradient of dynamic constraints
        for k = 0:N-1
            % Extract variables
            x_k = Y(3*k + 1);
            y_k = Y(3*k + 2);
            theta_k = Y(3*k + 3);
            
            % Indices in Y
            idx_x_k = 3*k + 1;
            idx_y_k = 3*k + 2;
            idx_theta_k = 3*k + 3;
            idx_x_k1 = 3*(k+1) + 1;
            idx_y_k1 = 3*(k+1) + 2;
            
            % Indices in ceq
            idx_ceq_x = 2*(k+1) + 1;
            idx_ceq_y = 2*(k+1) + 2;
            
            % Constraint for x[k+1]:
            % ceq_x = x[k+1] - x[k] - T*(w*cos(theta[k]) + y[k])
            
            % Partial derivatives of ceq_x
            gradceq(idx_x_k, idx_ceq_x) = -1;  % d/dx[k]
            gradceq(idx_y_k, idx_ceq_x) = -T;  % d/dy[k]
            gradceq(idx_theta_k, idx_ceq_x) = T * w * sin(theta_k);  % d/dtheta[k]
            gradceq(idx_x_k1, idx_ceq_x) = 1;  % d/dx[k+1]
            
            % Constraint for y[k+1]:
            % ceq_y = y[k+1] - y[k] - T*w*sin(theta[k])
            
            % Partial derivatives of ceq_y
            gradceq(idx_y_k, idx_ceq_y) = -1;  % d/dy[k]
            gradceq(idx_theta_k, idx_ceq_y) = -T * w * cos(theta_k);  % d/dtheta[k]
            gradceq(idx_y_k1, idx_ceq_y) = 1;  % d/dy[k+1]
        end
            %fprintf('Size of Jacobian: %d x %d\n', size(gradceq));
    end
end
