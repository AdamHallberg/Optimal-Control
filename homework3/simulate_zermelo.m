function [x, y] = simulate_zermelo(theta, Xi, T, w, N)
    % Forward simulate the Zermelo problem given control sequence
    x = zeros(N+1, 1);
    y = zeros(N+1, 1);
    
    % Initial condition
    x(1) = Xi(1);
    y(1) = Xi(2);
    
    % Euler integration
    for k = 1:N
        x(k+1) = x(k) + T * (w * cos(theta(k)) + y(k));
        y(k+1) = y(k) + T * w * sin(theta(k));
    end
end
