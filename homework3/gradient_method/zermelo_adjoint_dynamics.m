function dlambda_dt = zermelo_adjoint_dynamics(t, lambda, x, y, theta, t_span, w)
    x_t = interp1(t_span, x, t, 'linear', 'extrap');
    y_t = interp1(t_span, y, t, 'linear', 'extrap');
    theta_t = interp1(t_span, theta, t, 'linear', 'extrap');
    
    lambda_x = lambda(1);
    lambda_y = lambda(2);
    
    dlambda_x_dt = 0;
    dlambda_y_dt = -lambda_x;
    
    dlambda_dt = [dlambda_x_dt; dlambda_y_dt];
end
