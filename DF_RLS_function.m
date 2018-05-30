function [A,P] = DF_RLS_function(A,r_r,r_t,P)
    lambda = 0.9999; % Forgetting factor
    K = P*r_t;
    P = (1/lambda)*(P-((K*K')/(lambda+K'*r_t)));
    A = A-r_r*r_t'*P;
end

