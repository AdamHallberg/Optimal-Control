function [J, grad_J] = zermelo_objective_penalty(theta, Xi, T, w, N)

    [x, y] = simulate_zermelo(theta, Xi, T, w, N);
    
    J = -x(end);
    
    if nargout > 1
        lambda_x = zeros(N+1, 1);
        lambda_y = zeros(N+1, 1);
        
        lambda_x(N+1) = -1;
        lambda_y(N+1) = 0;
        
        for k = N:-1:1
            lambda_x(k) = lambda_x(k+1);
            lambda_y(k) = lambda_y(k+1) + T * lambda_x(k+1);
        end
        
        grad_J = zeros(N+1, 1);
        for k = 1:N
            grad_J(k) = -lambda_x(k) * w * sin(theta(k)) + lambda_y(k) * w * cos(theta(k));
        end
        grad_J(N+1) = 0;  % No control at final time
    end
end
