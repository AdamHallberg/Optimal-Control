function [F0, gradF0] = zermeloCostFun(Y, N)

    x_tf = Y(3*N + 1);
    
    F0 = -x_tf;
    
    if nargout > 1
        gradF0 = zeros(size(Y));
        gradF0(3*N + 1) = -1;
    end
end
