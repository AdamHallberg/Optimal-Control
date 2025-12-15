function [J, grad_J] = zermelo_objective_penalty(theta, Xi, T, w, N)
    % Objective function: minimize -x(tf)
    % Constraints are automatically satisfied by forward simulation
    
    [x, y] = simulate_zermelo(theta, Xi, T, w, N);
    
    % Cost: -x(tf)
    J = -x(end);
    
    % Compute gradient using adjoint method
    if nargout > 1
        % Adjoint equations (backward in time)
        lambda_x = zeros(N+1, 1);
        lambda_y = zeros(N+1, 1);
        
        % Terminal condition
        lambda_x(N+1) = -1;  % dJ/dx(tf) = -1
        lambda_y(N+1) = 0;   % dJ/dy(tf) = 0
        
        % Backward pass
        for k = N:-1:1
            % Adjoint equations (discrete-time)
            lambda_x(k) = lambda_x(k+1);
            lambda_y(k) = lambda_y(k+1) + T * lambda_x(k+1);
        end
        
        % Gradient w.r.t. control
        grad_J = zeros(N+1, 1);
        for k = 1:N
            % dH/dtheta = -lambda_x * w * sin(theta) + lambda_y * w * cos(theta)
            grad_J(k) = -lambda_x(k) * w * sin(theta(k)) + lambda_y(k) * w * cos(theta(k));
        end
        grad_J(N+1) = 0;  % No control at final time
    end
end
