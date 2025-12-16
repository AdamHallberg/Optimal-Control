function dXdt = zermelo_state_dynamics(t, X, theta, t_span, w)
    % State dynamics for Zermelo problem
    % X = [x; y]
    % 
    % Inputs:
    %   t      - current time
    %   X      - current state [x; y]
    %   theta  - control trajectory (vector for all time points)
    %   t_span - time points corresponding to theta
    %   w      - ship speed
    
    % Interpolate control at current time
    theta_t = interp1(t_span, theta, t, 'linear', 'extrap');
    
    % State equations from (1)
    % ẋ = w*cos(θ) + v(y) = w*cos(θ) + y
    % ẏ = w*sin(θ)
    
    x = X(1);
    y = X(2);
    
    dx_dt = w * cos(theta_t) + y;
    dy_dt = w * sin(theta_t);
    
    dXdt = [dx_dt; dy_dt];
end
