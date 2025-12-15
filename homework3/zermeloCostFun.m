function [F0, gradF0] = zermeloCostFun(Y, N)
    % Cost function for Zermelo problem: F0(Y) = -x(tf)
    % 
    % Inputs:
    %   Y - optimization variable vector [X[0]^T, theta[0], ..., X[N]^T, theta[N]]
    %   N - number of time steps
    %
    % Outputs:
    %   F0     - cost function value
    %   gradF0 - gradient of cost function
    
    % Extract x(tf) = x[N] which is at position 3*N+1 in Y
    x_tf = Y(3*N + 1);
    
    % Cost function: minimize -x(tf) to maximize x(tf)
    F0 = -x_tf;
    
    % Compute gradient
    if nargout > 1
        % Gradient is zero everywhere except at position 3*N+1
        gradF0 = zeros(size(Y));
        gradF0(3*N + 1) = -1;
    end
end
