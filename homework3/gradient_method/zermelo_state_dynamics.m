function dXdt = zermelo_state_dynamics(t, X, theta, t_span, w)
    theta_t = interp1(t_span, theta, t, 'linear', 'extrap');

    x = X(1);
    y = X(2);
    
    dx_dt = w * cos(theta_t) + y;
    dy_dt = w * sin(theta_t);
    
    dXdt = [dx_dt; dy_dt];
end
