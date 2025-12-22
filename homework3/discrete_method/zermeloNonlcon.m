function [c, ceq, gradc, gradceq] = zermeloNonlcon(Y, Xi, T, w, N)

    c = [];
    
    ceq = zeros(2*(N+1), 1);
    
    X0 = [Y(1); Y(2)];
    ceq(1:2) = X0 - Xi;
    
    % X[k+1] = F_bar(X[k], theta[k]) for k=0,...,N-1
    for k = 0:N-1
        x_k = Y(3*k + 1);
        y_k = Y(3*k + 2);
        theta_k = Y(3*k + 3);
        
        x_k1 = Y(3*(k+1) + 1);
        y_k1 = Y(3*(k+1) + 2);

        x_k1_pred = x_k + T * (w * cos(theta_k) + y_k);
        y_k1_pred = y_k + T * w * sin(theta_k);
        
        ceq(2*(k+1) + 1) = x_k1 - x_k1_pred;
        ceq(2*(k+1) + 2) = y_k1 - y_k1_pred;
    end
    
    if nargout > 2
        gradc = [];

        gradceq = zeros(length(Y), length(ceq));

        gradceq(1, 1) = 1;
        gradceq(2, 2) = 1;
        
        for k = 0:N-1
            x_k = Y(3*k + 1);
            y_k = Y(3*k + 2);
            theta_k = Y(3*k + 3);
            
            idx_x_k = 3*k + 1;
            idx_y_k = 3*k + 2;
            idx_theta_k = 3*k + 3;
            idx_x_k1 = 3*(k+1) + 1;
            idx_y_k1 = 3*(k+1) + 2;
            
            idx_ceq_x = 2*(k+1) + 1;
            idx_ceq_y = 2*(k+1) + 2;

            gradceq(idx_x_k, idx_ceq_x) = -1;
            gradceq(idx_y_k, idx_ceq_x) = -T;
            gradceq(idx_theta_k, idx_ceq_x) = T * w * sin(theta_k);
            gradceq(idx_x_k1, idx_ceq_x) = 1;

            gradceq(idx_y_k, idx_ceq_y) = -1;
            gradceq(idx_theta_k, idx_ceq_y) = -T * w * cos(theta_k);
            gradceq(idx_y_k1, idx_ceq_y) = 1;
        end
    end
end
