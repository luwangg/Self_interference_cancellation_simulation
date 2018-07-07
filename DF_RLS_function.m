function [a_RLS,P] = DF_RLS_function(a_RLS,y_RLS,x,P)
    lambda = 0.9999; % Forgetting factor
    K = P*x;
    P = (1/lambda)*(P-((K*K')/(lambda+K'*x)));
    a_RLS = a_RLS + conj(y_RLS)*P*x;
end