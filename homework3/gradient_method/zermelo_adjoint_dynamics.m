function dlambda_dt = zermelo_adjoint_dynamics(t, lambda, x, y, theta, t_span, w)
    % Adjoint dynamics for Zermelo problem
    % From Exercise 1.1(a)
    %
    % Inputs:
    %   t       - current time
    %   lambda  - current adjoint state [lambda_x; lambda_y]
    %   x       - state x trajectory (vector for all time points)
    %   y       - state y trajectory (vector for all time points)
    %   theta   - control trajectory (vector for all time points)
    %   t_span  - time points
    %   w       - ship speed
    
    % Interpolate at current time
    x_t = interp1(t_span, x, t, 'linear', 'extrap');
    y_t = interp1(t_span, y, t, 'linear', 'extrap');
    theta_t = interp1(t_span, theta, t, 'linear', 'extrap');
    
    lambda_x = lambda(1);
    lambda_y = lambda(2);
    
    % Adjoint equations from Exercise 1.1(a):
    % λ̇ = -∂H/∂X
    %
    % H = λ_x * (w*cos(θ) + y) + λ_y * w*sin(θ)
    %
    % ∂H/∂x = 0
    % ∂H/∂y = λ_x
    %
    % Therefore:
    % λ̇_x = -∂H/∂x = 0
    % λ̇_y = -∂H/∂y = -λ_x
    
    dlambda_x_dt = 0;
    dlambda_y_dt = -lambda_x;
    
    dlambda_dt = [dlambda_x_dt; dlambda_y_dt];
end
